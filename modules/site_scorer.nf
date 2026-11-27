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

// Optional prior from the protein alignments of the EviAnn run (bin/site_prior.py):
// start/stop probabilities at the CDS ends of miniprot alignments and/or EviAnn's
// best_unused_proteins are raised to params.site_prior_prob before UniAnn.
process SITE_PRIOR {
    tag "$id"

    input:
    tuple val(id), path(seq), path(sites)
    path run_dir

    output:
    tuple val(id), path('sites.prior.tsv'), emit: sites

    script:
    def g = params.genome.tokenize('/').last()
    def gffs = []
    if (params.site_prior in ['palign', 'both']) gffs << "--gff ${run_dir}/${g}.palign.fixed.gff --min-identity ${params.site_prior_min_identity}"
    if (params.site_prior in ['unused', 'both']) gffs << "--gff ${run_dir}/${g}.best_unused_proteins.gff"
    """
    site_prior.py ${seq} ${sites} ${gffs.join(' ')} --prob ${params.site_prior_prob} > sites.prior.tsv
    check_sites_tsv.py ${seq} sites.prior.tsv
    """
}
