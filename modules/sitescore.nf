// Site evaluator (sitescore submodule): donor/acceptor/start/stop probabilities.
// The model behind it (convmamba today) is chosen by params.sitescore_model.
// Output contract: docs/contracts.md §sites.tsv

process SITESCORE_TRAIN {
    publishDir "${params.outdir}/sitescore", mode: 'copy'

    input:
    path genome
    path annotation

    output:
    path 'model_dir', emit: model

    script:
    def init = params.sitescore_init ? "--init ${params.sitescore_init}" : ''
    def hp   = params.sitescore_hparams ? "--hparams '${params.sitescore_hparams}'" : ''
    """
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    sitescore train --model ${params.sitescore_model} --genome ${genome} --annotation ${annotation} \\
        --out model_dir ${init} ${hp}
    """
}

process SITESCORE_SCORE {
    tag "$id"

    input:
    tuple val(id), path(seq)
    path model_dir

    output:
    tuple val(id), path('sites.tsv'), emit: sites

    script:
    """
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    sitescore score --model-dir ${model_dir} --fasta ${seq} > sites.tsv
    """
}
