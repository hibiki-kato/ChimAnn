// One FASTA per segment: sequences >= params.min_seq_len are cut into
// params.segment_len pieces overlapping by params.segment_overlap
// (bin/segment_genome.py). Downstream PSAURON / site scoring / UniAnn work per
// segment (UniAnn needs ~0.6 GB RAM per Mb); UNIANN maps coordinates back and
// keeps each transcript from the segment that owns its midpoint.
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
