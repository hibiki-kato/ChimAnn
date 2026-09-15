// Alignment + filtering: the reference EviAnn driver, run inside its own
// directory so a later pass (INTEGRATE) can resume it with extra evidence.
process EVIANN {
    tag "$genome.name"
    publishDir "${params.outdir}/eviann", mode: 'copy', pattern: '*.gff'

    input:
    path genome

    output:
    path 'eviann_run',                      emit: run_dir
    path "${genome.name}.pseudo_label.gff", emit: gff

    script:
    def rna  = params.rnaseq   ? "-r ${params.rnaseq}"   : ''
    def prot = params.proteins ? "-p ${params.proteins}" : ''
    def src = params.eviann_src ? "export PATH=${params.eviann_src}:\$PATH" : ""
    """
    ${src}
    mkdir eviann_run && cd eviann_run
    eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args}
    cp ${genome.name}.pseudo_label.gff ..
    """
}
