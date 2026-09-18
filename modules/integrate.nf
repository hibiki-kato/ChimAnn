// Integration is EviAnn's job: rerun the driver on its own work directory with
// UniAnn's predictions passed as low-trust external CDS (-c + --untrusted-cds),
// so they go through EviAnn's splice-site filtering and final merge.
// Only predictions at loci with no evidence-based mRNA are passed (bin/novel_cds.py):
// EviAnn outputs external CDS only there anyway, and passing the rest makes the
// external alignments displace real protein evidence (measured on D. melanogaster:
// 631 correct CDS lost).
process INTEGRATE {
    tag "$genome.name"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'chimann.gff'

    input:
    path genome
    path run_dir
    path evidence_gff
    path ab_initio_gff

    output:
    path 'chimann.gff', emit: gff

    script:
    def rna  = params.rnaseq   ? "-r ${params.rnaseq}"   : ''
    def prot = params.proteins ? "-p ${params.proteins}" : ''
    def src = params.eviann_src ? "export PATH=${params.eviann_src}:\$PATH" : ""
    """
    ${src}
    novel_cds.py ${evidence_gff} ${ab_initio_gff} > novel.gff
    # -c only reruns merge and later, so the sorted BAMs (read by the completed
    # assembly stages only) are left out of the copy; hard links are unsafe here
    # because eviann.sh truncates existing files in place. Subdirectories
    # (the trained CNN splice model) are copied so the rerun reuses the model.
    mkdir run && find -L ${run_dir} -mindepth 1 -maxdepth 1 ! -name "*.bam" -exec cp -r -t run {} + && cd run
    eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args} \\
        -c ../novel.gff --untrusted-cds
    cp ${genome.name}.pseudo_label.gff ../chimann.gff
    """
}
