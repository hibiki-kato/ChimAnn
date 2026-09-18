// Integration is EviAnn's job: rerun the driver on its own work directory with
// UniAnn's predictions passed as low-trust external CDS (-c + --untrusted-cds),
// so they go through EviAnn's splice-site filtering and final merge. EviAnn
// (>= 2c1e1f9) keeps evidence-based CDS over external ones itself, so all
// predictions are passed; params.integrate_novel_only restores the older
// prefilter (bin/novel_cds.py) that kept only loci without evidence mRNA.
// k-best alternatives first go through bin/kbest_calibrate.py, which drops
// those whose score drop exceeds params.kbest_delta_fraction x the largest drop
// among alternatives confirmed by EviAnn's own alternative isoforms.
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
    def calib = (params.uniann_kbest as int) > 0 ?
        "kbest_calibrate.py ${evidence_gff} ${ab_initio_gff} --fraction ${params.kbest_delta_fraction} > calibrated.gff" :
        "cp ${ab_initio_gff} calibrated.gff"
    def novel = params.integrate_novel_only ? "novel_cds.py ${evidence_gff} calibrated.gff > novel.gff" : "cp calibrated.gff novel.gff"
    """
    ${src}
    ${calib}
    ${novel}
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
