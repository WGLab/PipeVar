nextflow.enable.dsl = 2

include { INPUT_PREPARATION } from '../../subworkflows/input_preparation/main'

process COMPATIBILITY_SMOKE {
    tag 'input-preparation-route'

    input:
    val route_settings

    output:
    path 'compatibility-smoke.txt', emit: result

    script:
    """
    printf '%s\n' '${route_settings.route}|${route_settings.batch_input}|${route_settings.effective_mode}' > compatibility-smoke.txt
    """
}

workflow {
    prepared = INPUT_PREPARATION()
    COMPATIBILITY_SMOKE(prepared.route_settings)
    COMPATIBILITY_SMOKE.out.result.view { result -> "SMOKE|${result.text.trim()}" }
}
