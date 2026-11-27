#!/usr/bin/env python3
"""Invoke CNN-Mamba training and oriented inference from an external checkout.

All model, loss and probability transformations remain in CNN-Mamba.
"""

import argparse
import hashlib
import json
import os
import shutil
import shlex
from pathlib import Path
import subprocess
import sys
import time


def read_sequence(fasta, chrom=None):
    records = []
    name, parts = None, []
    with Path(fasta).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(parts).upper()))
                name, parts = line[1:].split()[0], []
            elif name is not None:
                parts.append(line.strip())
        if name is not None:
            records.append((name, ''.join(parts).upper()))
    selected = [record for record in records if chrom is None or record[0] == chrom]
    if len(selected) != 1 or not selected[0][1]:
        raise ValueError('Select exactly one nonempty FASTA record with --chrom')
    return selected[0]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class CNNMamba:
    def __init__(self, repo, python=sys.executable):
        self.repo = Path(repo).resolve()
        self.python = str(python)
        self.env = dict(os.environ)
        self.env['PYTHONPATH'] = str(self.repo / 'src') + (
            os.pathsep + self.env['PYTHONPATH'] if self.env.get('PYTHONPATH') else '')
        self.commands = []

    def call(self, arguments, log=None):
        command = [self.python, *map(str, arguments)]
        started = time.time()
        if log is None:
            result = subprocess.run(command, env=self.env)
        else:
            with Path(log).open('w') as handle:
                result = subprocess.run(command, env=self.env, stdout=handle,
                                        stderr=subprocess.STDOUT)
        self.commands.append({'command': command, 'log': str(log) if log else None,
                              'exit_code': result.returncode,
                              'elapsed_seconds': time.time() - started})
        if result.returncode:
            if log is not None:
                print(Path(log).read_text()[-4000:], file=sys.stderr)
            raise subprocess.CalledProcessError(result.returncode, command)

    def module(self, module, arguments, log=None):
        self.call(['-m', 'cnn_mamba.' + module, *arguments], log)

    def train(self, arguments):
        self.module('train', arguments)

    def calibrate(self, arguments):
        self.module('calibrate_joint_frame', arguments)

    def model(self, fasta, annotation, outdir, checkpoint=None, calibration_json=None,
              validation_chrom=(), test_chrom=(), window_size=10000, stride=5000,
              train_args='', fit_temperature=False, include_test_in_train=False):
        if calibration_json and not checkpoint:
            raise ValueError('--calibration-json requires --checkpoint; fit temperature for a newly trained model')
        if checkpoint and fit_temperature:
            raise ValueError('--fit-temperature requires training data; supply --calibration-json for a reused model')
        if include_test_in_train and (checkpoint or fit_temperature):
            raise ValueError('--include-test-in-train applies to training without temperature fitting')
        if calibration_json and fit_temperature:
            raise ValueError('Select --calibration-json or --fit-temperature')
        outdir = Path(outdir).resolve()
        outdir.mkdir(parents=True, exist_ok=False)
        if checkpoint:
            shutil.copyfile(checkpoint, outdir/'best_model.pt')
        else:
            if not validation_chrom or not test_chrom:
                raise ValueError('Training requires explicit validation and test chromosomes')
            data = outdir.parent/'training_data'
            options = [Path(__file__).with_name('prepare_cnn_data.py'),
                       '--fasta', Path(fasta).resolve(), '--annotation', Path(annotation).resolve(),
                       '--outdir', data, '--window-size', window_size, '--stride', stride]
            for name in validation_chrom:
                options += ['--validation-chrom', name]
            for name in test_chrom:
                options += ['--test-chrom', name]
            self.call(options, outdir/'prepare.log')
            extras = shlex.split(train_args)
            # These alter data ownership or disable the required joint head.
            reserved = {'--data-root', '--output-dir', '--mode', '--window-size', '--stride',
                        '--all-splits-as-train', '--include-test-in-train',
                        '--joint-frame-aux-weight', '--phase-aux-weight', '--cds-aux-weight'}
            if any(arg.startswith('--') and any(option.startswith(arg.split('=')[0])
                                                for option in reserved) for arg in extras):
                raise ValueError('--train-args must preserve held-out splits and the joint four-class head')
            self.module('train', ['--data-root', data, '--mode', 'baseline',
                        '--output-dir', outdir, '--window-size', window_size, '--stride', stride,
                        '--epochs', 5, '--batch-size', 1, '--grad-accum', 8, '--workers', 0,
                        '--frame-dilation', 3, '--joint-frame-aux-weight', 0.1,
                        '--reference-name', 'EviAnn held-out pseudo-labels',
                        '--test-name', 'EviAnn held-out pseudo-labels',
                        # Test records then train too; validation still drives model selection.
                        *(['--include-test-in-train'] if include_test_in_train else []),
                        *extras], outdir/'train.log')
            shutil.copyfile(data/'metadata.json', outdir/'training-data.json')
            if fit_temperature:
                calibration = outdir.parent/'calibration'
                self.module('calibrate_joint_frame', [
                    '--checkpoint', outdir/'best_model.pt', '--data-root', data,
                    '--split', 'val', '--output-dir', calibration,
                    '--batch-size', 1, '--workers', 0], outdir/'calibrate.log')
                calibration_json = calibration/'temperature.json'
        if not (outdir/'best_model.pt').is_file():
            raise ValueError('CNN-Mamba produced no best_model.pt')
        if calibration_json:
            shutil.copyfile(calibration_json, outdir/'temperature.json')
        metadata = {'schema': 'chimann_cnn_model_v1', 'trained': not bool(checkpoint),
                    'checkpoint_sha256': sha256(outdir/'best_model.pt'),
                    'source_checkpoint': str(checkpoint) if checkpoint else None,
                    'joint_temperature': 'temperature.json' if calibration_json else None,
                    'temperature_scope': 'joint four-class logits only; site heads unchanged',
                    'label_source': ('Supplied checkpoint; training labels are not inferred'
                                     if checkpoint else
                                     'EviAnn pseudo-labels; held-out agreement is not gold accuracy'),
                    'commands': self.commands}
        (outdir/'model.json').write_text(json.dumps(metadata, indent=2)+'\n')
        return outdir

    def predict(self, fasta, checkpoint, outdir, chrom=None, strand='both',
                calibration_json=None, temperature=None, batch_size=1):
        if calibration_json is not None and temperature is not None:
            raise ValueError('Select calibration_json or temperature')
        fasta, checkpoint = Path(fasta).resolve(), Path(checkpoint).resolve()
        calibration_json = Path(calibration_json).resolve() if calibration_json else None
        name, sequence = read_sequence(fasta, chrom)
        outdir = Path(outdir).resolve()
        # Each request gets fresh output, avoiding upstream completion sentinels
        # silently reusing predictions made with a different model/temperature.
        outdir.mkdir(parents=True, exist_ok=False)
        manifest = {'schema': 'chimann_cnn_mamba_v1',
                    'input_fasta': str(fasta), 'input_fasta_sha256': sha256(fasta),
                    'checkpoint': str(checkpoint), 'checkpoint_sha256': sha256(checkpoint),
                    'calibration_json': str(calibration_json) if calibration_json else None,
                    'calibration_sha256': sha256(calibration_json) if calibration_json else None,
                    'repository': str(self.repo), 'python': self.python,
                    'sequence_name': name, 'sequence_length': len(sequence), 'jobs': []}
        manifest['repository_commit'] = subprocess.check_output(
            ['git', '-C', str(self.repo), 'rev-parse', 'HEAD'], text=True).strip()
        manifest['repository_local_changes'] = subprocess.check_output(
            ['git', '-C', str(self.repo), 'status', '--short'], text=True).splitlines()
        orientations = ('plus', 'minus') if strand == 'both' else (strand,)
        for orientation in orientations:
            job = outdir / orientation
            job.mkdir()
            bases = sequence if orientation == 'plus' else sequence.translate(
                str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN'))[::-1]
            oriented = job / 'sequence.fna'
            oriented.write_text('>' + name + '\n' + '\n'.join(
                bases[i:i + 80] for i in range(0, len(bases), 80)) + '\n')
            self.call(['-c', 'import pysam,sys; pysam.faidx(sys.argv[1])', oriented], job/'index.log')
            common = ['--checkpoint', checkpoint, '--fasta', oriented, '--chrom', name,
                      '--display-chrom', name, '--batch-size', str(batch_size)]
            print(f'{orientation}: scoring coding/frame and splice/start/stop', file=sys.stderr, flush=True)
            # Vector calibration files hold a per-class scale and bias; scalar ones one temperature.
            flag = '--vector-calibration-json' if calibration_json and 'vector_scale' in json.loads(
                calibration_json.read_text()) else '--calibration-json'
            joint_options = [flag, calibration_json] if calibration_json else (
                ['--temperature', str(temperature)] if temperature is not None else [])
            self.module('score_joint_frame_chromosome',
                        [*common, *joint_options, '--output-dir', job/'joint'], job/'joint.log')
            self.call([Path(__file__).with_name('score_cnn_sites.py'),
                       *common, '--output-dir', job/'sites'], job/'sites.log')
            self.module('export_joint_frame_tsv', [
                '--probabilities', job/'joint/joint_frame_probabilities.npy',
                '--argmax', job/'joint/joint_frame_argmax.npy',
                '--output', job/'coding.tsv', '--chrom', name, '--strand', '+'], job/'export-coding.log')
            self.module('export_uniann_scores', [
                '--score-dir', job/'sites', '--output', job/'sites.tsv',
                '--chrom', name, '--fasta', oriented, '--accession', name], job/'export-sites.log')
            metadata = json.loads((job/'joint/metadata.json').read_text())
            manifest['jobs'].append({
                'directory': str(job), 'original_strand': '+' if orientation == 'plus' else '-',
                'coordinate_system': 'one-based positions on sequence.fna; all TSV strand labels are +',
                'original_position_1based': 'p' if orientation == 'plus' else 'sequence_length-p+1',
                'fasta': str(oriented), 'coding': str(job/'coding.tsv'), 'sites': str(job/'sites.tsv'),
                'temperature': metadata['temperature'], 'temperature_scope': metadata['temperature_scope']})
        manifest['commands'] = self.commands
        manifest['completed'] = True
        result = outdir / 'manifest.json'
        result.write_text(json.dumps(manifest, indent=2) + '\n')
        return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True, help='CNN-Mamba checkout')
    parser.add_argument('--python', default=sys.executable, help='Python in the existing CNN-Mamba GPU environment')
    sub = parser.add_subparsers(dest='command', required=True)
    predict = sub.add_parser('predict', help='Export coding/frame and site probabilities on oriented strands')
    predict.add_argument('--fasta', type=Path, required=True)
    predict.add_argument('--checkpoint', type=Path, required=True)
    predict.add_argument('--outdir', type=Path, required=True, help='New output directory')
    predict.add_argument('--chrom')
    predict.add_argument('--strand', choices=('plus', 'minus', 'both'), default='both')
    predict.add_argument('--batch-size', type=int, default=1)
    group = predict.add_mutually_exclusive_group()
    group.add_argument('--calibration-json', type=Path)
    group.add_argument('--temperature', type=float)
    model = sub.add_parser('model', help='Train on EviAnn labels or bundle a supplied checkpoint')
    model.add_argument('--fasta', type=Path, required=True)
    model.add_argument('--annotation', type=Path, required=True)
    model.add_argument('--outdir', type=Path, required=True)
    model.add_argument('--checkpoint', type=Path)
    model.add_argument('--calibration-json', type=Path)
    model.add_argument('--validation-chrom', action='append', default=[])
    model.add_argument('--test-chrom', action='append', default=[])
    model.add_argument('--window-size', type=int, default=10000)
    model.add_argument('--stride', type=int, default=5000)
    model.add_argument('--train-args', default='')
    model.add_argument('--fit-temperature', action='store_true')
    model.add_argument('--include-test-in-train', action='store_true')
    for command in ('train', 'calibrate'):
        child = sub.add_parser(command, help='Forward existing CLI options after --; requires prepared data')
        child.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    client = CNNMamba(args.repo, args.python)
    try:
        if args.command == 'predict':
            output = client.predict(args.fasta, args.checkpoint, args.outdir,
                                    args.chrom, args.strand, args.calibration_json,
                                    args.temperature, args.batch_size)
            print(output)
        elif args.command == 'model':
            client.model(args.fasta, args.annotation, args.outdir, args.checkpoint,
                         args.calibration_json, args.validation_chrom, args.test_chrom,
                         args.window_size, args.stride, args.train_args, args.fit_temperature,
                         args.include_test_in_train)
        else:
            forwarded = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            getattr(client, args.command)(forwarded)
    except subprocess.CalledProcessError as error:
        return error.returncode
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
