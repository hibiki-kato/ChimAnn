// Model inference and all probability calibration stay in Bowen's implementation.
include { gitRevision; sourceDigest } from './eviann'

def cnnSourceDigest() {
    sourceDigest([file("${projectDir}/bin/run_cnn_mamba.py"), file("${projectDir}/bin/prepare_cnn_data.py"),
                  file("${projectDir}/bin/score_cnn_sites.py"), files("${params.cnn_mamba_dir}/src/cnn_mamba/*.py")])
}

def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
def chromOptions(option, values) {
    def names = values instanceof List ? values : (values ? values.toString().split(',').toList() : [])
    names.collect { "${option} ${quote(it)}" }.join(' ')
}

process CNN_MODEL {
    publishDir "${params.outdir}/cnn", mode: 'copy', pattern: 'model'
    input:
    path genome
    path annotation
    path checkpoints
    path calibrations
    val sources
    output:
    path 'model', emit: model
    script:
    def revision = gitRevision(params.cnn_mamba_dir)
    def weights = checkpoints ? "--checkpoint ${quote(checkpoints[0])}" : ''
    def temperature = calibrations ? "--calibration-json ${quote(calibrations[0])}" : ''
    def validation = chromOptions('--validation-chrom', params.cnn_validation_chroms)
    def test = chromOptions('--test-chrom', params.cnn_test_chroms)
    def fit = params.cnn_fit_temperature ? '--fit-temperature' : ''
    def no_test = params.cnn_include_test_in_train ? '--include-test-in-train' : ''
    """
    # cnn-mamba ${revision}
    # scoring sources ${sources}
    ${quote(params.cnn_python)} ${projectDir}/bin/run_cnn_mamba.py --repo ${quote(params.cnn_mamba_dir)} --python ${quote(params.cnn_python)} model \
        --fasta ${quote(genome)} --annotation ${quote(annotation)} --outdir model \
        ${weights} ${temperature} ${validation} ${test} ${fit} ${no_test} \
        --window-size ${params.cnn_window_size} --stride ${params.cnn_stride} --train-args ${quote(params.cnn_train_args)}
    """
}

process CNN_SCORE {
    tag "$id"
    publishDir { "${params.outdir}/cnn/scores/${id}" }, mode: 'copy', pattern: 'cnn_scores'
    input:
    tuple val(id), path(seq)
    path model
    val sources
    output:
    tuple val(id), path('cnn_scores'), emit: scores
    script:
    def revision = gitRevision(params.cnn_mamba_dir)
    """
    # cnn-mamba ${revision}
    # scoring sources ${sources}
    temperature=()
    if [ -f ${model}/temperature.json ]; then temperature=(--calibration-json ${model}/temperature.json); fi
    ${quote(params.cnn_python)} ${projectDir}/bin/run_cnn_mamba.py --repo ${quote(params.cnn_mamba_dir)} --python ${quote(params.cnn_python)} predict \
        --fasta ${quote(seq)} --checkpoint ${model}/best_model.pt --outdir cnn_scores --strand both \
        --batch-size ${params.cnn_batch_size} "\${temperature[@]}"
    """
}
