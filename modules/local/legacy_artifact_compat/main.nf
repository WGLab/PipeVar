// Materialize an nf-core caller output under PipeVar's public filename,
// decompressing bgzipped artifacts when the public contract requires it.
process MATERIALIZE_PUBLIC_ARTIFACT {
    tag "${meta.id}:${source_tool}"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14e1d96665f934a98e569fc5a6fa237f98d3753eee2b6f60d0aea8ff9d44f406/data' :
        'community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba' }"

    input:
    tuple val(meta), path(caller_output), val(source_tool), val(public_filename)

    output:
    tuple val(meta), path("${public_filename}"), val(source_tool), emit: artifact

    script:
    def materialize = caller_output.name.endsWith('.gz') \
    ? "bgzip --decompress --stdout \"${caller_output}\" > \"${public_filename}\"" \
    : "cp \"${caller_output}\" \"${public_filename}\""
    """
    ${materialize}
    """

    stub:
    """
    touch "${public_filename}"
    """
}
