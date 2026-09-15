// Integration is EviAnn's job: rerun the driver on its own work directory with
// UniAnn's predictions passed as low-trust external CDS (-c + --untrusted-cds),
// so they go through EviAnn's splice-site filtering and final merge.
process INTEGRATE {
    tag "$genome.name"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'chimann.gff'

    input:
    path genome
    path run_dir
    path ab_initio_gff

    output:
    path 'chimann.gff', emit: gff

    script:
    def rna  = params.rnaseq   ? "-r ${params.rnaseq}"   : ''
    def prot = params.proteins ? "-p ${params.proteins}" : ''
    def src = params.eviann_src ? "export PATH=${params.eviann_src}:\$PATH" : ""
    """
    ${src}
    cp -r ${run_dir}/ run && cd run
    eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args} \\
        -c ../${ab_initio_gff} --untrusted-cds
    cp ${genome.name}.pseudo_label.gff ../chimann.gff
    """
}
