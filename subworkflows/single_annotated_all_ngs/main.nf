include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Manta } from '../../modules/manta/'
include { scramble } from '../../modules/scramble/'
include { scramble_ref_prep } from '../../modules/scramble_ref_prep/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { CNVnator } from '../../modules/cnvnator/'
include { merge_shortread_sv_callers } from '../../modules/merge_shortread_sv_callers/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
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

// Single sample: imported ANNOVAR SNV + short-read SV/CNV/(optional) mito analysis.
workflow SINGLE_ANNOTATED_ALL_NGS {
	take:
	bam
	annovar_txt
	snv_vcf
	annovar_sv_vcf
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
	validator_script = Channel.value(file("${projectDir}/scripts/validate_preannotated_annovar_pair.py"))
	validate_preannotated_annovar_pair(
		out_prefix.combine(annovar_txt).combine(snv_vcf).map { prefix, annovar_txt_file, annovar_vcf_file ->
			tuple(prefix, annovar_txt_file, annovar_vcf_file)
		},
		validator_script
	)

	hpo = phenotype
	if ( phenotype_format == 'clinical_note' || phenotype_format == 'yes' ) {
		phenotagger(phenotype, out_prefix)
		hpo = phenotagger.out
	}
	Phen2gene(hpo, out_prefix)
	RankVar(validate_preannotated_annovar_pair.out[0], Phen2gene.out, hpo, out_prefix, gnomad, gq, ad, rankvar_filter)
	rankscore_result = Rankscore_analysis(validate_preannotated_annovar_pair.out[0], Phen2gene.out, out_prefix, gnomad, rankscore_filter, rankscore_softwares, gq, phen2gene_top_n)

	ANNOVAR_SV(annovar_sv_vcf, out_prefix, Phen2gene.out, "null", "preannotated")
	SURVIVOR(ANNOVAR_SV.out, out_prefix)
	PhenoSV(SURVIVOR.out, out_prefix, hpo)
	ExpansionHunter(bam, out_prefix, eh_ref_fa, eh_variant_catalog)
	eh_filter(out_prefix, ExpansionHunter.out.json)
	ngs_prio(out_prefix, RankVar.out, rankscore_result.out.rankscore, rankscore_result.out.clinvar, PhenoSV.out, ANNOVAR_SV.out, validate_preannotated_annovar_pair.out[1], hpo, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	if ( mito_ref_fa != null ) {
		mutect2_ref_fa = mito_ref_fa.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
			tuple(fa_file, fai_file, dict_file)
		}
		mito_prep_mutect2(bam, out_prefix, mito_ref_fa, mito_contig)
		mito_mutect2(mito_prep_mutect2.out, out_prefix, mutect2_ref_fa, mito_contig)
		mito_annotation(mito_mutect2.out[0], out_prefix)
		mito_prio(mito_annotation.out[0], out_prefix)
	}

	emit:
	prio_vcf = ngs_prio.out[0]
	prio_gene_vcf = ngs_prio.out[1]
	mito_vcf = mito_ref_fa != null ? mito_mutect2.out[0] : Channel.empty()
	mito_annotated_tsv = mito_ref_fa != null ? mito_annotation.out[0] : Channel.empty()
	mito_prioritized_tsv = mito_ref_fa != null ? mito_prio.out : Channel.empty()
}
