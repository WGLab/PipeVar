include { multi_denovo_snv_vcf_filter } from '../../modules/multi_denovo_snv_vcf_filter/'

// Cohort-level raw/called SNV filtering before phenotype-specific annotation.
workflow DENOVO_SNV_VCF_FILTER_CORE {
	take:
	snv_records
	pedigree_csv
	role_column
	family_column
	vcf_sample_column
	exclude_contigs

	main:
	cohort_input = snv_records.collect().map { rows ->
		def ordered = rows.sort { left, right -> left[0].toString() <=> right[0].toString() }
		tuple(ordered.collect { it[0] }, ordered.collect { it[1] })
	}
	filtered = multi_denovo_snv_vcf_filter(
		cohort_input,
		pedigree_csv,
		role_column,
		family_column,
		vcf_sample_column,
		exclude_contigs
	)
	vcf_files = filtered.vcf.flatten().collect()
	filtered_records = filtered.bindings
		.splitCsv(header: true, sep: '\t')
		.combine(vcf_files)
		.map { binding, all_vcf ->
			def vcf_file = all_vcf.find { it.name == binding.vcf_path }
			if (vcf_file == null) {
				throw new IllegalStateException("de novo raw SNV binding references missing output for '${binding.out_prefix}'")
			}
			tuple(binding.out_prefix, vcf_file)
		}

	emit:
	records = filtered_records
	summary = filtered.summary
	bindings = filtered.bindings
}
