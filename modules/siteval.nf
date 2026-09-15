// Site evaluator (siteval submodule): donor/acceptor/start/stop probabilities.
// The model behind it (SSM today) is chosen by params.siteval_model.
// Output contract: docs/contracts.md §sites.tsv

process SITEVAL_TRAIN {
    publishDir "${params.outdir}/siteval", mode: 'copy'

    input:
    path genome
    path annotation

    output:
    path 'model_dir', emit: model

    script:
    def init = params.siteval_init ? "--init ${params.siteval_init}" : ''
    def hp   = params.siteval_hparams ? "--hparams '${params.siteval_hparams}'" : ''
    """
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    siteval train --model ${params.siteval_model} --genome ${genome} --annotation ${annotation} \\
        --out model_dir ${init} ${hp}
    """
}

process SITEVAL_SCORE {
    tag "$id"

    input:
    tuple val(id), path(seq)
    path model_dir

    output:
    tuple val(id), path('sites.tsv'), emit: sites

    script:
    """
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    siteval score --model-dir ${model_dir} --fasta ${seq} > sites.tsv
    """
}
