#!/usr/bin/env python3
"""Keep UniAnn models and add distinct, complete AltAnn alternatives.

Both GFFs must use the same genomic coordinates.
Every AltAnn rank-zero CDS chain must occur in UniAnn, including its strand and CDS phases.
Additional UniAnn models are retained.
AltAnn rank-zero rows are omitted, and rank-positive models are added only when their CDS chain is absent from all retained models.
Alternatives whose locus has no retained AltAnn reference are omitted because the reference failed an ownership or short-CDS filter.
AltAnn's escaped sequence identifiers are decoded once to match the original FASTA identifiers.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from urllib.parse import unquote


TRANSCRIPTS = {'transcript', 'mRNA'}


@dataclass
class Feature:
    line: str
    fields: list
    attrs: dict
    parents: tuple
    label: str

    @property
    def identifier(self):
        return self.attrs.get('ID')


class Annotation:
    def __init__(self, path, decode_seqid=False):
        self.rows = []
        self.ids = {}
        self.models = {}
        with Path(path).open() as handle:
            for number, line in enumerate(handle, 1):
                line = line.rstrip('\r\n')
                if not line or line.startswith('#'):
                    continue
                label = f'{path}:{number}'
                fields = line.split('\t')
                if len(fields) != 9:
                    raise ValueError(f'{label}: expected nine GFF columns')
                if decode_seqid:
                    fields[0] = unquote(fields[0])
                    line = '\t'.join(fields)
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise ValueError(f'{label}: invalid GFF coordinates') from error
                if start < 1 or end < start or fields[6] not in ('+', '-'):
                    raise ValueError(f'{label}: invalid coordinates or strand')
                if fields[2] == 'CDS' and fields[7] not in ('0', '1', '2'):
                    raise ValueError(f'{label}: CDS phase must be 0, 1, or 2')
                attrs = {}
                for item in fields[8].split(';'):
                    if not item or item == '.':
                        continue
                    if '=' not in item:
                        raise ValueError(f'{label}: malformed GFF attribute')
                    key, value = item.split('=', 1)
                    if not key or not value or key in attrs:
                        raise ValueError(f'{label}: empty or duplicate GFF attribute')
                    attrs[key] = value
                parents = tuple(attrs.get('Parent', '').split(',')) if 'Parent' in attrs else ()
                row = Feature(line, fields, attrs, parents, label)
                if row.identifier:
                    if row.identifier in self.ids:
                        raise ValueError(f'{label}: duplicate ID {row.identifier}')
                    self.ids[row.identifier] = row
                if fields[2] in TRANSCRIPTS:
                    if not row.identifier:
                        raise ValueError(f'{label}: transcript requires an ID')
                    self.models[row.identifier] = row
                elif fields[2] != 'gene' and not parents:
                    raise ValueError(f'{label}: child feature requires a Parent')
                self.rows.append(row)
        self._owners = {}
        self.chains = {identifier: [] for identifier in self.models}
        for row in self.rows:
            for parent in row.parents:
                if parent not in self.ids:
                    raise ValueError(f'{row.label}: missing Parent {parent}')
                ancestor = self.ids[parent]
                if row.fields[2] in TRANSCRIPTS and ancestor.fields[2] != 'gene':
                    raise ValueError(f'{row.label}: transcript Parent must be a gene')
                if (row.fields[0], row.fields[6]) != (ancestor.fields[0], ancestor.fields[6]):
                    raise ValueError(f'{row.label}: Parent has a different sequence or strand')
                if not (int(ancestor.fields[3]) <= int(row.fields[3]) <= int(row.fields[4])
                        <= int(ancestor.fields[4])):
                    raise ValueError(f'{row.label}: child lies outside its Parent')
            owners = self.owners(row)
            if row.fields[2] not in TRANSCRIPTS | {'gene'} and not owners:
                raise ValueError(f'{row.label}: child does not belong to a transcript')
            if row.fields[2] == 'CDS':
                for identifier in owners:
                    self.chains[identifier].append((int(row.fields[3]), int(row.fields[4]),
                                                    row.fields[7]))
        for identifier, chain in self.chains.items():
            row = self.models[identifier]
            chain.sort(reverse=row.fields[6] == '-')
            if len(chain) != len(set(chain)):
                raise ValueError(f'{row.label}: duplicate CDS interval')

    def owners(self, row, visiting=None):
        key = id(row)
        if key in self._owners:
            return self._owners[key]
        visiting = set() if visiting is None else visiting
        if key in visiting:
            raise ValueError(f'{row.label}: cyclic Parent relationship')
        visiting.add(key)
        parents = set()
        for parent in row.parents:
            if parent not in self.ids:
                raise ValueError(f'{row.label}: missing Parent {parent}')
            parents.update(self.owners(self.ids[parent], visiting))
        visiting.remove(key)
        owners = {row.identifier} if row.fields[2] in TRANSCRIPTS else parents
        self._owners[key] = owners
        return owners

    def signature(self, identifier):
        model = self.models[identifier]
        return model.fields[0], model.fields[6], tuple(self.chains[identifier])


def merge(uniann_path, altann_path):
    """Validate both annotations before returning any output rows."""
    original, alternatives = Annotation(uniann_path), Annotation(altann_path, decode_seqid=True)
    seen = {original.signature(identifier) for identifier in original.models
            if original.chains[identifier]}
    ranks = {}
    for identifier, model in alternatives.models.items():
        raw_rank = model.attrs.get('kbest_rank', '')
        if not raw_rank.isascii() or not raw_rank.isdigit():
            raise ValueError(f'{model.label}: AltAnn transcript requires a nonnegative kbest_rank')
        ranks[identifier] = int(raw_rank)
        chain = alternatives.chains[identifier]
        if not chain:
            raise ValueError(f'{model.label}: AltAnn transcript has no CDS')
        length = 0
        for start, end, phase in chain:
            if int(phase) != (-length) % 3:
                raise ValueError(f'{model.label}: inconsistent CDS phases')
            length += end - start + 1
        if length % 3:
            raise ValueError(f'{model.label}: AltAnn transcript is not a complete coding model')
        if ranks[identifier] == 0 and alternatives.signature(identifier) not in seen:
            raise ValueError(f'{model.label}: AltAnn reference CDS chain is absent from UniAnn')
    reference_loci = {(row.fields[0], row.fields[6], row.attrs['kbest_locus'])
                      for identifier, row in alternatives.models.items()
                      if ranks[identifier] == 0 and 'kbest_locus' in row.attrs}
    for identifier, row in alternatives.models.items():
        if ranks[identifier] > 0 and 'kbest_locus' not in row.attrs:
            raise ValueError(f'{row.label}: AltAnn alternative requires kbest_locus')
    for row in alternatives.rows:
        owners = alternatives.owners(row)
        if len(owners) > 1:
            raise ValueError(f'{row.label}: AltAnn feature belongs to multiple transcript blocks')
        if owners and 'kbest_rank' in row.attrs:
            rank = row.attrs['kbest_rank']
            if not rank.isascii() or not rank.isdigit() or int(rank) != ranks[next(iter(owners))]:
                raise ValueError(f'{row.label}: mixed ranks in an AltAnn transcript block')
    selected = set()
    for identifier in alternatives.models:
        model = alternatives.models[identifier]
        signature = alternatives.signature(identifier)
        locus = (model.fields[0], model.fields[6], model.attrs.get('kbest_locus'))
        if ranks[identifier] > 0 and locus in reference_loci and signature not in seen:
            selected.add(identifier)
            seen.add(signature)
    kept = {id(row) for row in alternatives.rows if alternatives.owners(row) & selected}
    for row in alternatives.rows:
        if id(row) not in kept:
            continue
        pending = list(row.parents)
        while pending:
            parent = alternatives.ids[pending.pop()]
            if id(parent) not in kept:
                kept.add(id(parent))
                pending.extend(parent.parents)
    output = [row for row in alternatives.rows if id(row) in kept]
    ids = set(original.ids)
    for row in output:
        if row.identifier:
            if row.identifier in ids:
                raise ValueError(f'{row.label}: output ID collides with UniAnn ID {row.identifier}')
            ids.add(row.identifier)
    return ['##gff-version 3', *[row.line for row in original.rows], *[row.line for row in output]]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('uniann_gff', help='Original UniAnn GFF in genomic coordinates')
    parser.add_argument('altann_gff', help='AltAnn GFF in the same coordinates with kbest_rank attributes')
    args = parser.parse_args(argv)
    try:
        rows = merge(args.uniann_gff, args.altann_gff)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    sys.stdout.write('\n'.join(rows) + '\n')


if __name__ == '__main__':
    main()
