// PSAURON coding-potential scores per sequence and strand (-a required by UniAnn).
// The '-' run gets the reverse-complemented sequence from REVCOMP.
process PSAURON {
    tag "$id $strand"

    input:
    tuple val(id), val(strand), path(seq)

    output:
    tuple val(id), val(strand), path('psauron_score.csv'), emit: csv

    script:
    // attempt 2 (after a CUDA OOM) falls back to CPU
    def dev = task.attempt > 1 ? 'CUDA_VISIBLE_DEVICES=""' : ''
    """
    export OMP_NUM_THREADS=${task.cpus}
    ${dev} psauron -i ${seq} -a
    """
}
