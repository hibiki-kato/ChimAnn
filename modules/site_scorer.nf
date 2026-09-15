// Site scorer stage: donor/acceptor/start/stop probabilities for UniAnn.
// Any tool can fill this slot; the pipeline only runs two command templates
// (params.site_train_cmd / params.site_score_cmd) and reads back sites.tsv.
// Contract: README "Bring your own site scorer" and docs/contracts.md.
// Defaults call the bundled sitescore submodule.

def fill(String tpl, Map v) { v.inject(tpl) { s, k, val -> s.replace('{' + k + '}', "${val}") } }

process SITE_TRAIN {
    publishDir "${params.outdir}/site_scorer", mode: 'copy'

    input:
    path genome
    path annotation
    val  gpu_barrier      // unused; forces ordering after other GPU tasks

    output:
    path 'model_dir', emit: model

    script:
    def init = params.sitescore_init ? "--init ${params.sitescore_init}" : ''
    def hp   = params.sitescore_hparams ? "--hparams '${params.sitescore_hparams}'" : ''
    def tpl  = params.site_train_cmd ?:
        "sitescore train --model ${params.sitescore_model} --genome {genome} --annotation {annotation} --out {model_dir} ${init} ${hp}"
    """
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    ${fill(tpl, [genome: genome, annotation: annotation, model_dir: 'model_dir'])}
    test -d model_dir
    """
}

process SITE_SCORE {
    tag "$id"

    input:
    tuple val(id), path(seq)
    path model_dir

    output:
    tuple val(id), path('sites.tsv'), emit: sites

    script:
    def tpl = params.site_score_cmd ?: "sitescore score --model-dir {model_dir} --fasta {fasta}"
    """
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    ${fill(tpl, [model_dir: model_dir, fasta: seq])} > sites.tsv
    check_sites_tsv.py ${seq} sites.tsv
    """
}
