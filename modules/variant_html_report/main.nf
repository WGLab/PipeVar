// Render self-contained final variant HTML reports from sample-keyed records.
// Small-variant/SV evidence only.
process variant_html_report_no_repeat {
    container = 'beoungl/docker_test:html_report'

    input:
    tuple val(out_prefix), path(prio_vcf), path(prio_gene_report)

    output:
    tuple val(out_prefix), path("${out_prefix}.variant_html_report.html")

    script:
	"""
	python3 /opt/pipevar/html_report/render_final_ngs_report.py \
	    --sample "$out_prefix" \
	    --prio-vcf "$prio_vcf" \
	    --prio-gene-vcf "$prio_gene_report" \
	    --output-html "${out_prefix}.variant_html_report.html" \
	    --title "PipeVar final variant report" \
	    --mito-min-af "${params.mito_gui_min_af}" \
	    --mito-min-apogee2 "${params.mito_gui_min_apogee2}" \
	    --mito-min-mitotip "${params.mito_gui_min_mitotip}"
	"""
}

// Add repeat-expansion evidence to the core report.
process variant_html_report_with_repeat {
    container = 'beoungl/docker_test:html_report'

    input:
    tuple val(out_prefix), path(prio_vcf), path(prio_gene_report), path(repeat_tsv)

    output:
    tuple val(out_prefix), path("${out_prefix}.variant_html_report.html")

    script:
	"""
	python3 /opt/pipevar/html_report/render_final_ngs_report.py \
	    --sample "$out_prefix" \
	    --prio-vcf "$prio_vcf" \
	    --prio-gene-vcf "$prio_gene_report" \
	    --repeat-tsv "$repeat_tsv" \
	    --output-html "${out_prefix}.variant_html_report.html" \
	    --title "PipeVar final variant report" \
	    --mito-min-af "${params.mito_gui_min_af}" \
	    --mito-min-apogee2 "${params.mito_gui_min_apogee2}" \
	    --mito-min-mitotip "${params.mito_gui_min_mitotip}"
	"""
}

// Add repeat-expansion and mitochondrial evidence to the core report.
process variant_html_report_with_repeat_and_mito {
    container = 'beoungl/docker_test:html_report'

    input:
    tuple val(out_prefix), path(prio_vcf), path(prio_gene_report), path(repeat_tsv), path(mito_tsv)

    output:
    tuple val(out_prefix), path("${out_prefix}.variant_html_report.html")

    script:
	"""
	python3 /opt/pipevar/html_report/render_final_ngs_report.py \
	    --sample "$out_prefix" \
	    --prio-vcf "$prio_vcf" \
	    --prio-gene-vcf "$prio_gene_report" \
	    --repeat-tsv "$repeat_tsv" \
	    --mito-tsv "$mito_tsv" \
	    --output-html "${out_prefix}.variant_html_report.html" \
	    --title "PipeVar final variant report" \
	    --mito-min-af "${params.mito_gui_min_af}" \
	    --mito-min-apogee2 "${params.mito_gui_min_apogee2}" \
	    --mito-min-mitotip "${params.mito_gui_min_mitotip}"
	"""
}
