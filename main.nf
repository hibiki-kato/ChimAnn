#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EVIANN; eviannSourceDigest } from './modules/eviann'
include { SPLIT_GENOME } from './modules/split_genome'
include { CNN_MODEL; CNN_SCORE; cnnSourceDigest } from './modules/cnn_mamba'
include { UNIANN; uniannSourceDigest } from './modules/uniann'
include { ALTANN; altannSourceDigest } from './modules/altann'
include { INTEGRATE    } from './modules/integrate'

workflow {
    def legacy = params.keySet().findAll { it.startsWith('psauron_') || it.startsWith('site_') || it.startsWith('sitescore_') || it == 'uniann_args' }
    if (legacy)
        error "Legacy scoring options ${legacy.join(', ')} are unsupported; use the CNN parameters in example/params.yaml"
    if (params.eviann_reuse_annotation && !params.eviann_resume_dir)
        error 'eviann_reuse_annotation requires eviann_resume_dir'
    if (!params.cnn_checkpoint && (!params.cnn_validation_chroms || !params.cnn_test_chroms))
        error 'Training requires cnn_validation_chroms and cnn_test_chroms (disjoint FASTA record lists)'
    genome = Channel.fromPath(params.genome, checkIfExists: true).first()
    eviann_sources = Channel.value(eviannSourceDigest())
    EVIANN(genome, eviann_sources)
    checkpoints = Channel.value(params.cnn_checkpoint ? [file(params.cnn_checkpoint, checkIfExists: true)] : [])
    calibrations = Channel.value(params.cnn_calibration_json ? [file(params.cnn_calibration_json, checkIfExists: true)] : [])
    cnn_sources = Channel.value(cnnSourceDigest())
    CNN_MODEL(genome, EVIANN.out.gff, checkpoints, calibrations, cnn_sources)
    SPLIT_GENOME(genome)
    seqs = SPLIT_GENOME.out.seqs.flatten().map { f -> tuple(f.baseName, f) }
    CNN_SCORE(seqs, CNN_MODEL.out.model, cnn_sources)
    UNIANN(CNN_SCORE.out.scores, Channel.value(uniannSourceDigest()))
    def k = params.uniann_kbest != null ? params.uniann_kbest as int : params.altann_kbest as int
    if (k > 0) {
        ALTANN(UNIANN.out.jobs, Channel.value(altannSourceDigest()))
        predictions = ALTANN.out.gff
    } else {
        predictions = UNIANN.out.jobs.map { id, job, gff -> tuple(id, gff) }
    }
    ab_initio_gff = predictions.map { id, gff -> gff }.collectFile(name: 'ab_initio.gff', sort: true)
    INTEGRATE(genome, EVIANN.out.run_dir, EVIANN.out.gff, ab_initio_gff, eviann_sources)
}
