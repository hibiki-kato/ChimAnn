// Integration is EviAnn's job: rerun the driver on its own work directory with
// AltAnn best paths and alternatives passed as external CDS (-c).
// EviAnn applies its integration rules and final merge. EviAnn
// (>= 2c1e1f9) keeps evidence-based CDS over external ones itself, so all
// predictions are passed; params.integrate_novel_only restores the older
// prefilter (bin/novel_cds.py) that kept only loci without evidence mRNA.
// Setting kbest_delta_fraction enables the optional evidence-calibrated score cutoff.
include { gitRevision } from './eviann'

process INTEGRATE {
    tag "$genome.name"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'chimann.gff'

    input:
    path genome
    path run_dir
    path evidence_gff
    path ab_initio_gff
    val sources

    output:
    path 'chimann.gff', emit: gff

    script:
    def rna  = params.rnaseq   ? "-r ${params.rnaseq}"   : ''
    def prot = params.proteins ? "-p ${params.proteins}" : ''
    def src = params.eviann_src ? "export PATH=${params.eviann_src}:\$PATH" : ""
    def rev = gitRevision(params.eviann_src)
    def calib = params.kbest_delta_fraction != null ?
        "kbest_calibrate.py ${evidence_gff} ${ab_initio_gff} --fraction ${params.kbest_delta_fraction} > calibrated.gff" :
        "cp ${ab_initio_gff} calibrated.gff"
    def novel = params.integrate_novel_only ? "novel_cds.py ${evidence_gff} calibrated.gff > novel.gff" : "cp calibrated.gff novel.gff"
    """
    # eviann ${rev}
    # eviann sources ${sources}
    ${src}
    ${calib}
    ${novel}
    # -c only reruns merge and later, so the sorted BAMs (read by the completed
    # assembly stages only) are left out of the copy; hard links are unsafe here
    # because eviann.sh truncates existing files in place. Subdirectories needed
    # by the completed run are copied as well.
    mkdir run && find -L ${run_dir} -mindepth 1 -maxdepth 1 ! -name "*.bam" -exec cp -rL -t run {} + && cd run
    # A failed driver must not publish the copied first-pass annotation as new output.
    rm -f ${genome.name}.pseudo_label.gff
    eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args} \\
        -c ../novel.gff ${params.integrate_args}
    n=\$(awk -F'\\t' '\$3=="mRNA"' ${genome.name}.pseudo_label.gff | wc -l)
    [ "\$n" -gt 0 ] || { echo "EviAnn integration produced no mRNA; see run/*.err" >&2; exit 1; }
    cp ${genome.name}.pseudo_label.gff ../chimann.gff
    """
}
