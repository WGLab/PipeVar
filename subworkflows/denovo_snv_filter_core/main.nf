include { multi_denovo_snv_filter } from '../../modules/multi_denovo_snv_filter/'

// Cohort-level de novo SNV filtering with explicit out_prefix/file bindings.
workflow DENOVO_SNV_FILTER_CORE {
	take:
	annovar_records
	pedigree_csv
	role_column
	family_column
	vcf_sample_column
	exclude_contigs

	main:
	cohort_input = annovar_records.collect().map { rows ->
		def ordered = rows.sort { left, right -> left[0].toString() <=> right[0].toString() }
		tuple(
			ordered.collect { it[0] },
			ordered.collect { it[1] },
			ordered.collect { it[2] }
		)
	}
	filtered = multi_denovo_snv_filter(
		cohort_input,
		pedigree_csv,
		role_column,
		family_column,
		vcf_sample_column,
		exclude_contigs
	)
	txt_files = filtered.txt.flatten().collect()
	vcf_files = filtered.vcf.flatten().collect()
	filtered_records = filtered.bindings
		.splitCsv(header: true, sep: '\t')
		.combine(txt_files)
		.combine(vcf_files)
		.map { binding, all_txt, all_vcf ->
			def txt_file = all_txt.find { it.name == binding.txt_path }
			def vcf_file = all_vcf.find { it.name == binding.vcf_path }
			if (txt_file == null || vcf_file == null) {
				throw new IllegalStateException("de novo SNV binding references missing output for '${binding.out_prefix}'")
			}
			tuple(binding.out_prefix, txt_file, vcf_file)
	}

	emit:
	records = filtered_records
	summary = filtered.summary
	bindings = filtered.bindings
}
