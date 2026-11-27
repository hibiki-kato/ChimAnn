// Decode exactly the CNN probabilities; retain standardized scores for AltAnn.
include { gitRevision; sourceDigest } from './eviann'

def uniannSourceDigest() {
    def decoder = params.ab_initio_decoder == 'altann'
        ? file("${params.altann_dir}/build/altann-core") : file("${params.uniann_dir}/bin/uniann")
    sourceDigest([file("${projectDir}/bin/run_uniann.py"),
        file("${params.uniann_dir}/scripts/uniann.sh"), file("${params.uniann_dir}/scripts/convert_coding_frame.pl"),
        decoder])
}

process UNIANN {
    tag "$id"
    publishDir { "${params.outdir}/uniann/${id}" }, mode: 'copy'
    // Measured peak of the native decoder with UniAnn PR #4: 0.205 GB per Mb. Without
    // that PR it needs about 0.31 GB per Mb, so long chromosomes would not fit.
    memory {
        def length = (id =~ /__(\d+)-(\d+)-\d+$/)[0].with { (it[2] as long) - (it[1] as long) }
        def perMb = params.ab_initio_decoder == 'uniann' ? 0.21 : 0.02
        "${Math.ceil(length / 1e6 * perMb + 4) as long} GB"
    }
    input:
    tuple val(id), path(scores)
    val sources
    output:
    tuple val(id), path('uniann_job'), path('uniann.gff3'), emit: jobs
    script:
    def revision = gitRevision(params.uniann_dir)
    def decoder = params.ab_initio_decoder == 'altann' ? "--altann-dir ${params.altann_dir}"
        : (params.uniann_trace ? '' : '--no-trace')
    """
    # uniann ${revision}
    # decoder ${params.ab_initio_decoder} ${params.ab_initio_decoder == 'altann' ? gitRevision(params.altann_dir) : ''}
    # decoder sources ${sources}
    export OMP_NUM_THREADS=${task.cpus}
    python3 ${projectDir}/bin/run_uniann.py --input ${scores} --uniann-dir ${params.uniann_dir} ${decoder} \
        --outdir uniann_job --output uniann.gff3 --mult ${params.score_multiplier} \
        --overlap ${params.segment_overlap} --margin ${params.segment_margin}
    """
}
