// One FASTA per sequence; downstream tools (PSAURON, UniAnn) work per sequence.
// Sequences shorter than params.min_seq_len skip the ab initio track (they are
// still annotated by EviAnn): each one would cost a full psauron model load.
process SPLIT_GENOME {
    input:
    path genome

    output:
    path 'seqs/*.fa', emit: seqs

    script:
    """
    mkdir seqs
    awk '/^>/{f="seqs/"substr(\$1,2)".fa"} {print > f}' ${genome}
    for f in seqs/*.fa; do
        len=\$(grep -v '^>' \$f | tr -d '\\n' | wc -c)
        [ "\$len" -lt ${params.min_seq_len} ] && rm \$f
    done
    ls seqs/*.fa >/dev/null
    """
}
