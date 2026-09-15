#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ChimAnn: EviAnn (alignment + filtering) -> PSAURON + siteval (UniAnn inputs)
//          -> UniAnn (ab initio) -> EviAnn again, taking UniAnn CDS as low-trust evidence

include { EVIANN        } from './modules/eviann'
include { SPLIT_GENOME  } from './modules/split_genome'
include { PSAURON       } from './modules/psauron'
include { SITEVAL_TRAIN } from './modules/siteval'
include { SITEVAL_SCORE } from './modules/siteval'
include { REVCOMP       } from './modules/uniann'
include { SPLIT_SITES   } from './modules/uniann'
include { UNIANN        } from './modules/uniann'
include { INTEGRATE     } from './modules/integrate'

workflow {
    genome = Channel.fromPath(params.genome, checkIfExists: true)

    // 1. Evidence-based annotation (EviAnn)
    EVIANN(genome)
    evidence_gff = EVIANN.out.gff

    // 2. Per-sequence scatter: UniAnn and PSAURON take one sequence at a time
    SPLIT_GENOME(genome)
    seqs = SPLIT_GENOME.out.seqs.flatten().map { f -> tuple(f.baseName, f) }

    // 3a. Both strands: UniAnn decodes + only, so "-" = reverse complement
    REVCOMP(seqs)
    stranded = seqs.map { id, f -> tuple(id, '+', f) }.mix(REVCOMP.out.seq)   // (id, strand, fa)

    // 3b. PSAURON coding-potential emissions per strand
    PSAURON(stranded)

    // 3c. Site evaluator: train on EviAnn preliminary annotation, score both strands, split
    SITEVAL_TRAIN(genome, evidence_gff)
    SITEVAL_SCORE(seqs, SITEVAL_TRAIN.out.model.collect())   // value channel: one model, every sequence
    SPLIT_SITES(seqs.join(SITEVAL_SCORE.out.sites))
    sites = SPLIT_SITES.out.plus.mix(SPLIT_SITES.out.minus)                 // (id, strand, tsv)

    // 4. Ab initio prediction (UniAnn), one run per sequence and strand
    UNIANN(stranded.join(PSAURON.out.csv, by: [0, 1]).join(sites, by: [0, 1]))
    ab_initio_gff = UNIANN.out.gff.map { it[2] }.collectFile(name: 'uniann.gff', sort: true)

    // 5. Integration by EviAnn: UniAnn CDS enter as low-trust external evidence
    INTEGRATE(genome, EVIANN.out.run_dir, ab_initio_gff)
}
