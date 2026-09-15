// One FASTA per sequence; downstream tools (PSAURON, UniAnn) work per sequence.
process SPLIT_GENOME {
    input:
    path genome

    output:
    path 'seqs/*.fa', emit: seqs

    script:
    """
    mkdir seqs
    awk '/^>/{f="seqs/"substr(\$1,2)".fa"} {print > f}' ${genome}
    """
}
