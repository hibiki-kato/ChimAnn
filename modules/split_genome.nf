// One FASTA per sequence >= params.min_seq_len, kept whole by default. A nonzero
// params.segment_len cuts pieces overlapping by params.segment_overlap
// (bin/segment_genome.py), and ALTANN keeps each transcript from the segment
// that owns its midpoint.
// Shorter sequences skip the ab initio track (EviAnn still annotates them).
process SPLIT_GENOME {
    input:
    path genome

    output:
    path 'seqs/*.fa', emit: seqs

    script:
    """
    segment_genome.py segment ${genome} --size ${params.segment_len} --overlap ${params.segment_overlap} \\
        --min-len ${params.min_seq_len} --outdir seqs
    """
}
