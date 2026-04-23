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
	out_prefix
	ref_fa
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
	if ( phenotype_format == 'clinical_note' ) {
		phenotagger(phenotype, out_prefix)
		hpo = phenotagger.out
	}
	Phen2gene(hpo, out_prefix)
	RankVar(validate_preannotated_annovar_pair.out[0], Phen2gene.out, hpo, out_prefix, gnomad, gq, ad, rankvar_filter)
	rankscore_result = Rankscore_analysis(validate_preannotated_annovar_pair.out[0], Phen2gene.out, out_prefix, gnomad, rankscore_filter, rankscore_softwares, gq, phen2gene_top_n)

	Manta(bam, out_prefix, ref_fa)
	scramble_mode = params.scramble ? params.scramble.toString().trim().toLowerCase() : "no"
	scramble_vcf = null
	if ( scramble_mode == "yes" ) {
		scramble_ref_meta = ref_fa.map { ref_tuple -> tuple([id: 'reference'], ref_tuple[0], ref_tuple[1]) }
		scramble_ref_bundle = scramble_ref_prep(scramble_ref_meta)
		scramble_cluster_input = out_prefix.combine(bam).map { prefix, bam_tuple ->
			tuple([id: prefix], bam_tuple[0], bam_tuple[1])
		}
		scramble(scramble_cluster_input, scramble_ref_bundle.out.ref)
		scramble_vcf = scramble.out.vcf.map { meta, vcf -> vcf }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && scramble_mode == "yes" ) {
		normalize_shortread_alignment(bam, out_prefix, ref_fa)
		CNVnator(normalize_shortread_alignment.out, out_prefix, ref_fa, params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(scramble_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs, out_prefix)
		sv_vcf = merge_shortread_sv_callers.out
	}
	else if ( cnvnator_mode != "no" ) {
		normalize_shortread_alignment(bam, out_prefix, ref_fa)
		CNVnator(normalize_shortread_alignment.out, out_prefix, ref_fa, params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs, out_prefix)
		sv_vcf = merge_shortread_sv_callers.out
	}
	else if ( scramble_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(scramble_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs, out_prefix)
		sv_vcf = merge_shortread_sv_callers.out
	}
	else {
		sv_vcf = Manta.out
	}

	ANNOVAR_SV(sv_vcf, out_prefix, Phen2gene.out, "null")
	SURVIVOR(ANNOVAR_SV.out, out_prefix)
	PhenoSV(SURVIVOR.out, out_prefix, hpo)
	ExpansionHunter(bam, out_prefix, ref_fa)
	eh_filter(out_prefix, ExpansionHunter.out)
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
