// Alignment + filtering: the reference EviAnn driver, run inside its own
// directory so a later pass (INTEGRATE) can resume it with extra evidence.
// params.eviann_resume_dir seeds the run from a finished EviAnn directory
// (its non-BAM files are copied) so the read/protein alignment stages are
// skipped; merge.success is cleared so filtering and merging rerun with the
// current eviann.sh and splice model.
// eviann_reuse_annotation keeps the completed annotation and skips this first-pass driver.
// INTEGRATE still reruns merging with the current driver and external CDS.
// Record the tool revision in the task script. CNN/decoder stages additionally
// pass content fingerprints as explicit inputs for -resume invalidation.
def gitRevision(dir) {
    if (!dir) return 'PATH'
    def p = ['git', '-C', "${dir}", 'describe', '--always', '--dirty'].execute()
    p.waitFor()
    p.exitValue() == 0 ? p.text.trim() : 'unknown'
}

// Absolute and nested Python helpers are not hashed as Nextflow bundled scripts.
def sourceDigest(paths) {
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    paths.flatten().sort { it.toString() }.each { digest.update(new File(it.toString()).bytes) }
    digest.digest().encodeHex().toString()
}

def eviannSourceDigest() {
    if (!params.eviann_src) return 'PATH'
    sourceDigest([files("${params.eviann_src}/*.sh"), files("${params.eviann_src}/*.pl"),
                  files("${params.eviann_src}/my_maker*")])
}

process EVIANN {
    tag "$genome.name"
    publishDir "${params.outdir}/eviann", mode: 'copy', pattern: '*.gff'

    input:
    path genome
    val sources

    output:
    path 'eviann_run',                      emit: run_dir
    path "${genome.name}.pseudo_label.gff", emit: gff

    script:
    def rna  = params.rnaseq   ? "-r ${params.rnaseq}"   : ''
    def prot = params.proteins ? "-p ${params.proteins}" : ''
    def src = params.eviann_src ? "export PATH=${params.eviann_src}:\$PATH" : ""
    def rev = gitRevision(params.eviann_src)
    def reset = params.eviann_reuse_annotation ? '' : ' && rm -f eviann_run/merge.success'
    def seed = params.eviann_resume_dir ?
        "find -L ${params.eviann_resume_dir} -mindepth 1 -maxdepth 1 ! -name '*.bam' -exec cp -rL -t eviann_run {} +${reset}" : ""
    def run = params.eviann_reuse_annotation ? '' :
        "eviann.sh -t ${task.cpus} -g ../${genome} ${rna} ${prot} ${params.eviann_args}"
    """
    # eviann ${rev}
    # eviann sources ${sources}
    ${src}
    mkdir eviann_run
    ${seed}
    cd eviann_run
    ${run}
    # eviann.sh can exit 0 after an internal tool failure; refuse an empty annotation
    n=\$(awk -F'\\t' '\$3=="mRNA"' ${genome.name}.pseudo_label.gff | wc -l)
    [ "\$n" -gt 0 ] || { echo "EviAnn produced no mRNA; see eviann_run/*.err" >&2; exit 1; }
    cp ${genome.name}.pseudo_label.gff ..
    """
}
