// Alignment + filtering: the reference EviAnn driver, run inside its own
// directory so a later pass (INTEGRATE) can resume it with extra evidence.
process EVIANN {
    tag "$genome.name"
    publishDir "${params.outdir}/eviann", mode: 'copy', pattern: 'eviann_run/*.pseudo_label.gff'

    input:
    path genome

    output:
    path 'eviann_run',                                 emit: run_dir
    path "eviann_run/${genome.name}.pseudo_label.gff", emit: gff

    script:
    def rna  = params.rnaseq   ? "-r ${params.rnaseq}"   : ''
    def prot = params.proteins ? "-p ${params.proteins}" : ''
    """
    mkdir eviann_run && cd eviann_run
    eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args}
    """
}
