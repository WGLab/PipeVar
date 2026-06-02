include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Manta } from '../../modules/manta/'
include { xtea } from '../../modules/xtea/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { CNVnator } from '../../modules/cnvnator/'
include { truvari_shortread_sv_merge } from '../../modules/truvari_shortread_sv_merge/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { eh_filter } from '../../modules/eh_filter/'
include { ngs_prio } from '../../modules/ngs_prio/'
include { validate_preannotated_annovar_pair } from '../../modules/validate_preannotated_annovar_pair/'
include { mito_prep_mutect2 } from '../../modules/mito_prep_mutect2/'
include { mito_mutect2 } from '../../modules/mito_mutect2/'
include { mito_annotation } from '../../modules/mito_annotation/'
include { mito_prio } from '../../modules/mito_prio/'
include { variant_html_report; variant_html_report_with_mito } from '../../modules/variant_html_report/'

// Single sample: imported ANNOVAR SNV + called short-read SV/CNV/(optional) mito analysis.
workflow SINGLE_ANNOTATED_SNV_CALLED_SV_NGS {
	take:
	bam
	annovar_txt
	snv_vcf
	out_prefix
	ref_fa
	eh_ref_fa
	eh_variant_catalog
	phenotype
	phenotype_format
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
	validated_annovar = validate_preannotated_annovar_pair(
		out_prefix.combine(annovar_txt).combine(snv_vcf).map { prefix, annovar_txt_file, annovar_vcf_file ->
			tuple(prefix, annovar_txt_file, annovar_vcf_file)
		}
	)
	validated_annovar_txt = validated_annovar.map { prefix, validated_txt, validated_vcf -> validated_txt }
	validated_annovar_vcf = validated_annovar.map { prefix, validated_txt, validated_vcf -> validated_vcf }

	hpo = phenotype
	if ( phenotype_format == 'clinical_note' || phenotype_format == 'yes' ) {
		phenotagger(phenotype, out_prefix)
		hpo = phenotagger.out
	}
	Phen2gene(hpo, out_prefix)
	RankVar(validated_annovar_txt, Phen2gene.out, hpo, out_prefix, gnomad, gq, ad, rankvar_filter)
	rankscore_result = Rankscore_analysis(validated_annovar_txt, Phen2gene.out, out_prefix, gnomad, rankscore_filter, rankscore_softwares, gq, phen2gene_top_n)

	Manta(bam, out_prefix, ref_fa)
	xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : "no"
	xtea_vcf = null
	if ( xtea_mode == "yes" ) {
		xtea_input = out_prefix.combine(bam).map { prefix, bam_tuple ->
			tuple([id: prefix], bam_tuple[0], bam_tuple[1])
		}
		xtea(xtea_input, ref_fa)
		xtea_vcf = xtea.out.vcf.map { meta, vcf -> vcf }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && xtea_mode == "yes" ) {
		normalize_shortread_alignment(bam, out_prefix, ref_fa)
		CNVnator(normalize_shortread_alignment.out, out_prefix, ref_fa, params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(xtea_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs, ref_fa, out_prefix)
		sv_vcf = truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( cnvnator_mode != "no" ) {
		normalize_shortread_alignment(bam, out_prefix, ref_fa)
		CNVnator(normalize_shortread_alignment.out, out_prefix, ref_fa, params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs, ref_fa, out_prefix)
		sv_vcf = truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( xtea_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(xtea_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs, ref_fa, out_prefix)
		sv_vcf = truvari_shortread_sv_merge.out.merged_vcf
	}
	else {
		sv_vcf = Manta.out
	}

	ANNOVAR_SV(sv_vcf, out_prefix, Phen2gene.out, "null", "called")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out, out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream, out_prefix)
	PhenoSV(SURVIVOR.out, out_prefix, hpo)
	ExpansionHunter(bam, out_prefix, eh_ref_fa, eh_variant_catalog)
	eh_filter(out_prefix, ExpansionHunter.out)
	ngs_prio(out_prefix, RankVar.out, rankscore_result.out.rankscore, rankscore_result.out.clinvar, PhenoSV.out, annovar_sv_for_downstream, validated_annovar_vcf, hpo, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	if ( mito_ref_fa != null ) {
		mutect2_ref_fa = mito_ref_fa.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
			tuple(fa_file, fai_file, dict_file)
		}
		mito_prep_mutect2(bam, out_prefix, mito_ref_fa, mito_contig)
		mito_mutect2(mito_prep_mutect2.out, out_prefix, mutect2_ref_fa, mito_contig)
		mito_annotation(mito_mutect2.out[0], out_prefix)
		mito_prio(mito_annotation.out[0], out_prefix)
		variant_html_report_with_mito(out_prefix, ngs_prio.out[0], ngs_prio.out[1], eh_filter.out, mito_prio.out)
	}
	else {
		variant_html_report(out_prefix, ngs_prio.out[0], ngs_prio.out[1], eh_filter.out)
	}

	emit:
	prio_vcf = ngs_prio.out[0]
	prio_gene_vcf = ngs_prio.out[1]
	mito_vcf = mito_ref_fa != null ? mito_mutect2.out[0] : Channel.empty()
	mito_annotated_tsv = mito_ref_fa != null ? mito_annotation.out[0] : Channel.empty()
	mito_prioritized_tsv = mito_ref_fa != null ? mito_prio.out : Channel.empty()
}
