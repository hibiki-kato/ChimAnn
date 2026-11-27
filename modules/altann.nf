// Reuse UniAnn's emissions, transitions and native traceback without recalibration.
include { gitRevision; sourceDigest } from './eviann'

def altannSourceDigest() {
    sourceDigest([file("${projectDir}/bin/run_altann.py"),
        files("${params.altann_dir}/altann/*.py"), file("${params.altann_dir}/build/altann-core")])
}

process ALTANN {
    tag "$id"
    publishDir { "${params.outdir}/altann/${id}" }, mode: 'copy', pattern: '*.gff*'
    input:
    tuple val(id), path(job), path(uniann_gff, stageAs: 'uniann_best.gff3')
    val sources
    output:
    tuple val(id), path("${id}.ab_initio.gff"), emit: gff
    tuple val(id), path("${id}.altann.gff3"), path("${id}.altann.gff3.tsv"), path("${id}.altann.gff3.json"), emit: alternatives
    script:
    def revision = gitRevision(params.altann_dir)
    def k = params.uniann_kbest != null ? params.uniann_kbest as int : params.altann_kbest as int
    """
    # altann ${revision}
    # decoder sources ${sources}
    python3 ${projectDir}/bin/run_altann.py --input ${job} --altann-dir ${params.altann_dir} \
        --output ${id}.altann.gff3 --k ${k} --flank ${params.altann_flank} --threads ${task.cpus} \
        --overlap ${params.segment_overlap} --margin ${params.segment_margin}
    merge_altann.py ${uniann_gff} ${id}.altann.gff3 > ${id}.ab_initio.gff
    """
}
