// Decode the CNN-Mamba probabilities with AltAnn alone (UniAnn removed 2026-11-27):
// temperature-free inputs from `altann prepare`, the global best path from exact
// K=1 search (--exact-baseline) and local K-best alternatives.
include { gitRevision; sourceDigest } from './eviann'

def altannSourceDigest() {
    sourceDigest([file("${projectDir}/bin/run_altann.py"), file("${projectDir}/bin/merge_altann.py"),
        files("${params.altann_dir}/altann/*.py"), file("${params.altann_dir}/build/altann-core")])
}

process ALTANN {
    tag "$id"
    publishDir { "${params.outdir}/altann/${id}" }, mode: 'copy', pattern: '*.gff*'
    input:
    tuple val(id), path(scores)
    val sources
    output:
    tuple val(id), path("${id}.ab_initio.gff"), emit: gff
    tuple val(id), path("${id}.altann.gff3"), path("${id}.altann.gff3.tsv"), path("${id}.altann.gff3.json"), emit: alternatives
    script:
    def revision = gitRevision(params.altann_dir)
    def exact = params.altann_metadata_states ? '--metadata-states' : ''
    """
    # altann ${revision}
    # decoder sources ${sources}
    python3 ${projectDir}/bin/run_altann.py --input ${scores} --altann-dir ${params.altann_dir} \
        --output ${id}.altann.gff3 --ab-initio ${id}.ab_initio.gff --k ${params.altann_kbest} \
        --flank ${params.altann_flank} --threads ${task.cpus} --mult ${params.score_multiplier} \
        --overlap ${params.segment_overlap} --margin ${params.segment_margin} ${exact}
    """
}
