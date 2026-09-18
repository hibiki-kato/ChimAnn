// Alignment + filtering: the reference EviAnn driver, run inside its own
// directory so a later pass (INTEGRATE) can resume it with extra evidence.
// params.eviann_resume_dir seeds the run from a finished EviAnn directory
// (its non-BAM files are copied) so the read/protein alignment stages are
// skipped; merge.success is cleared so filtering and merging rerun with the
// current eviann.sh and splice model.
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
    def seed = params.eviann_resume_dir ?
        "find -L ${params.eviann_resume_dir} -mindepth 1 -maxdepth 1 ! -name '*.bam' -exec cp -r -t eviann_run {} + && rm -f eviann_run/merge.success" : ""
    """
    ${src}
    mkdir eviann_run
    ${seed}
    cd eviann_run
    eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args}
    # eviann.sh can exit 0 after an internal tool failure; refuse an empty annotation
    n=\$(awk -F'\\t' '\$3=="mRNA"' ${genome.name}.pseudo_label.gff | wc -l)
    [ "\$n" -gt 0 ] || { echo "EviAnn produced no mRNA; see eviann_run/*.err" >&2; exit 1; }
    cp ${genome.name}.pseudo_label.gff ..
    """
}
