// Ab initio gene prediction with UniAnn on one segment and one strand; the output
// is mapped back to genome coordinates (segment_genome.py unsegment).
// UniAnn decodes the + strand only, so "-" runs on the reverse complement
// (REVCOMP) with the sites table remapped, and the GFF is flipped back here.
process REVCOMP {
    tag "$id"

    input:
    tuple val(id), path(seq)

    output:
    tuple val(id), val('-'), path("${id}.rc.fa"), emit: seq

    script:
    """
    strand_tools.py revcomp ${seq} > ${id}.rc.fa
    """
}

process SPLIT_SITES {
    tag "$id"

    input:
    tuple val(id), path(seq), path(sites)

    output:
    tuple val(id), val('+'), path('sites.plus.tsv'), emit: plus
    tuple val(id), val('-'), path('sites.rc.tsv'),   emit: minus

    script:
    """
    awk -F'\\t' 'NR==1 || \$3=="+"' ${sites} > sites.plus.tsv
    strand_tools.py sites_rc ${seq} ${sites} > sites.rc.tsv
    """
}

process UNIANN {
    tag "$id $strand"
    publishDir "${params.outdir}/uniann", mode: 'copy'

    input:
    tuple val(id), val(strand), path(seq), path(psauron_csv), path(sites)

    output:
    tuple val(id), val(strand), path("${id}.${strand == '+' ? 'plus' : 'minus'}.uniann.gff"), emit: gff

    script:
    def out = "${id}.${strand == '+' ? 'plus' : 'minus'}.uniann.gff"
    def back = strand == '+' ? "cat ${seq}.uniann.gff" : "strand_tools.py flip_gff ${seq} ${seq}.uniann.gff | gffread"
    back += " | segment_genome.py unsegment /dev/stdin --overlap ${params.segment_overlap}"
    """
    export OMP_NUM_THREADS=${task.cpus}
    # uniann.sh scales scores by the best donor and dies unless max(donor prob) * e > 1;
    # sequences with no confident donor (e.g. mitochondria) get an empty annotation.
    if awk -F'\t' 'NR>1 && \$4=="donor" && \$NF*2.718281828 > 1 {found=1} END {exit !found}' ${sites}; then
        ${params.uniann_dir}/bin/uniann.sh -f ${seq} -p ${psauron_csv} -s ${sites} ${params.uniann_args}
        ${back} > ${out}
    else
        echo "no donor with prob > 1/e in ${sites}; skipping UniAnn on ${id} ${strand}" >&2
        printf '##gff-version 3\n' > ${out}
    fi
    """
}
