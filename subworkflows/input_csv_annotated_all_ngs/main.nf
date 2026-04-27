include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_scramble } from '../../modules/multi_scramble/'
include { scramble_ref_prep } from '../../modules/scramble_ref_prep/'
include { multi_normalize_shortread_alignment } from '../../modules/multi_normalize_shortread_alignment/'
include { multi_cnvnator } from '../../modules/multi_cnvnator/'
include { multi_merge_shortread_sv_callers } from '../../modules/multi_merge_shortread_sv_callers/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_ngs_prio } from '../../modules/multi_ngs_prio/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { validate_preannotated_annovar_pair } from '../../modules/validate_preannotated_annovar_pair/'
include { multi_mito_prep_mutect2 } from '../../modules/multi_mito_prep_mutect2/'
include { multi_mito_mutect2 } from '../../modules/multi_mito_mutect2/'
include { multi_mito_annotation } from '../../modules/multi_mito_annotation/'
include { multi_mito_prio } from '../../modules/multi_mito_prio/'

// CSV batch: imported ANNOVAR SNV + short-read SV/CNV/(optional) mito analysis.
workflow INPUT_CSV_ANNOTATED_ALL_NGS {
	take:
	input_annotated_ngs
	input_age
	ref_fa
	eh_ref_fa
	eh_variant_catalog
	rankscore_filter
	rankscore_softwares
	phen2gene_top_n
	gnomad
	gq
	ad
	rankvar_filter
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet
	mito_ref_fa
	mito_contig

	main:
	validator_script = Channel.value(file("${projectDir}/scripts/validate_preannotated_annovar_pair.py"))
	validate_input = input_annotated_ngs.map { out_prefix, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
		tuple(out_prefix, annovar_txt, annovar_vcf)
	}
	validated_annovar = validate_preannotated_annovar_pair(validate_input, validator_script)

	clinical_note_input = input_annotated_ngs
		.filter { out_prefix, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> phenotype_format == 'clinical_note' }
		.map { out_prefix, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }
	hpo_input = input_annotated_ngs
		.filter { out_prefix, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> phenotype_format == 'hpo' }
		.map { out_prefix, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }
	phenotagger_result = multi_phenotagger(clinical_note_input)
	hpo_paths = phenotagger_result.mix(hpo_input)
	phen2gene_result = multi_phen2gene(hpo_paths)

	validated_annovar_txt = validated_annovar.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_txt) }
	join_annovar_phen2gene = validated_annovar_txt.join(phen2gene_result)
	join_annovar_hpo = join_annovar_phen2gene.join(hpo_paths)
	rankscore_result = multi_rankscore(join_annovar_phen2gene, gnomad, rankscore_filter, rankscore_softwares, gq, phen2gene_top_n)
	rankvar_result = multi_rankvar(join_annovar_hpo, gnomad, gq, ad, rankvar_filter)

	input_bam_with_bam = input_annotated_ngs.map { out_prefix, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
		tuple(out_prefix, bam_file, bai_file)
	}
	multi_eh_result = multi_expansionhunter(input_bam_with_bam, eh_ref_fa, eh_variant_catalog)
	multi_eh_filter(multi_eh_result.json)
	manta_result = multi_manta(input_bam_with_bam, ref_fa)

	scramble_mode = params.scramble ? params.scramble.toString().trim().toLowerCase() : "no"
	scramble_vcf = null
	if ( scramble_mode == "yes" ) {
		scramble_ref_meta = ref_fa.map { ref_tuple -> tuple([id: 'reference'], ref_tuple[0], ref_tuple[1]) }
		scramble_ref_bundle = scramble_ref_prep(scramble_ref_meta)
		scramble_cluster_input = input_bam_with_bam.map { out_prefix, bam_file, index_file ->
			tuple([id: out_prefix], bam_file, index_file)
		}
		multi_scramble(scramble_cluster_input, scramble_ref_bundle.out.ref)
		scramble_vcf = multi_scramble.out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && scramble_mode == "yes" ) {
		normalized_bam = multi_normalize_shortread_alignment(input_bam_with_bam, ref_fa)
		multi_cnvnator(normalized_bam, ref_fa, params.cnvnator_bin_size)
		merged_sv_input = manta_result.join(multi_cnvnator.out.vcf).join(scramble_vcf).map { out_prefix, manta_vcf, cnvnator_vcf, scramble_vcf_file ->
			tuple(out_prefix, [manta_vcf, cnvnator_vcf, scramble_vcf_file])
		}
		sv_result = multi_merge_shortread_sv_callers(merged_sv_input)
	}
	else if ( cnvnator_mode != "no" ) {
		normalized_bam = multi_normalize_shortread_alignment(input_bam_with_bam, ref_fa)
		multi_cnvnator(normalized_bam, ref_fa, params.cnvnator_bin_size)
		merged_sv_input = manta_result.join(multi_cnvnator.out.vcf).map { out_prefix, manta_vcf, cnvnator_vcf ->
			tuple(out_prefix, [manta_vcf, cnvnator_vcf])
		}
		sv_result = multi_merge_shortread_sv_callers(merged_sv_input)
	}
	else if ( scramble_mode == "yes" ) {
		merged_sv_input = manta_result.join(scramble_vcf).map { out_prefix, manta_vcf, scramble_vcf_file ->
			tuple(out_prefix, [manta_vcf, scramble_vcf_file])
		}
		sv_result = multi_merge_shortread_sv_callers(merged_sv_input)
	}
	else {
		sv_result = manta_result
	}

	sv_result_annovar = sv_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file ->
		tuple(out_prefix, vcf_file, phen2gene_file, "null")
	}
	annovar_sv_result = multi_annovar_sv(sv_result_annovar)
	survivor_result = multi_survivor(annovar_sv_result)
	phenosv_input = survivor_result.join(hpo_paths)
	phenosv_result = multi_phenosv(phenosv_input)

	validated_annovar_vcf = validated_annovar.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_vcf) }
	phenosv_annovar_snv = phenosv_result.join(validated_annovar_vcf)
	sv_join = phenosv_annovar_snv.join(annovar_sv_result)
	rankscore_join = sv_join.join(rankscore_result)
	rankvar_join = rankscore_join.join(rankvar_result)
	hpo_with_age = hpo_paths.join(input_age).map { out_prefix, hpo_path, age_of_onset -> tuple(out_prefix, hpo_path, age_of_onset) }
	rankvar_join_hpo = rankvar_join.join(hpo_with_age)
	rankvar_join_hpo_ordered = rankvar_join_hpo.map { out_prefix, sv_pathogenic, snv_vcf_path, sv_vcf_path, snv_rankscore, snv_pathogenic, snv_rankvar, hpo_path, age_of_onset ->
		tuple(out_prefix, snv_rankvar, snv_rankscore, snv_pathogenic, sv_pathogenic, sv_vcf_path, snv_vcf_path, hpo_path, age_of_onset)
	}
	multi_ngs_prio(rankvar_join_hpo_ordered, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	if ( mito_ref_fa != null ) {
		mutect2_ref_fa = mito_ref_fa.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
			tuple(fa_file, fai_file, dict_file)
		}
		prepped = multi_mito_prep_mutect2(input_bam_with_bam, mito_ref_fa, mito_contig)
		mt_vcf = multi_mito_mutect2(prepped, mutect2_ref_fa, mito_contig)
		annotated = multi_mito_annotation(mt_vcf)
		multi_mito_prio(annotated)
	}

	emit:
	prio_vcf = multi_ngs_prio.out[0]
	prio_gene_vcf = multi_ngs_prio.out[1]
	mito_vcf = mito_ref_fa != null ? mt_vcf : Channel.empty()
	mito_annotated_tsv = mito_ref_fa != null ? annotated : Channel.empty()
	mito_prioritized_tsv = mito_ref_fa != null ? multi_mito_prio.out : Channel.empty()
}
