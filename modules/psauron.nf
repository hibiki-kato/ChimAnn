// PSAURON coding-potential scores for all six frames of each original sequence.
// run in chunks (bin/psauron_chunked.py). CPU only: psauron is faster on 4 CPU
// threads than on the GPU (measured 33 s vs 40 s per 4.5 Mb) and this keeps the
// GPU free for the site scorer; tasks overlap with EviAnn.
process PSAURON {
    tag "$id"

    input:
    tuple val(id), path(seq)

    output:
    tuple val(id), path('psauron_score.csv'), emit: csv

    script:
    """
    export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=${task.cpus}
    psauron_chunked.py ${seq} --chunk ${params.psauron_chunk} --out psauron_score.csv
    """
}
