

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Manta } from '../../modules/manta/'
include { scramble } from '../../modules/scramble/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { CNVnator } from '../../modules/cnvnator/'
include { merge_shortread_sv_callers } from '../../modules/merge_shortread_sv_callers/'
include { deepvariant } from '../../modules/deepvariant/'
include { haplotypecaller } from '../../modules/haplotypecaller/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { eh_filter } from '../../modules/eh_filter/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { ngs_prio } from '../../modules/ngs_prio/'


// Single sample: short-read full path (SNP + SV + STR) with DeepVariant and Manta.
workflow SINGLE_ALIGNMENT_ALL_NGS {
	take:
	bam
	out_prefix
	ref_fa
	note
	rankscore_filter
	rankscore_softwares
	phen2gene_top_n
	gnomad
	gq
	ad

	rankvar_filter
	is_note
	target
	caller_mode
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet

	main:

	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	Phen2gene(hpo,out_prefix)
	if ( target == "yes" ) {
		phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix,phen2gene_top_n)
		if ( caller_mode == "haplotypecaller" ) {
			haplotypecaller(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
		else {
			deepvariant(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
	}
	else {
		if ( caller_mode == "haplotypecaller" ) {
			haplotypecaller(bam,out_prefix,ref_fa,target)
		}
		else {
			deepvariant(bam,out_prefix,ref_fa,target)
		}
	}
	snp_vcf = (caller_mode == "haplotypecaller") ? haplotypecaller.out : deepvariant.out
	annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR(snp_vcf,out_prefix,annovar_bed)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	Manta(bam,out_prefix,ref_fa)
	scramble_mode = params.scramble ? params.scramble.toString().trim().toLowerCase() : "no"
	scramble_vcf = null
	if ( scramble_mode == "yes" ) {
		scramble_ref_meta = ref_fa.map { ref_tuple -> tuple([id: 'reference'], ref_tuple[0], ref_tuple[1]) }
		scramble_cluster_input = out_prefix.combine(bam).map { prefix, bam_tuple ->
			tuple([id: prefix], bam_tuple[0], bam_tuple[1])
		}
		scramble(scramble_cluster_input, scramble_ref_meta)
		scramble_vcf = scramble.out.vcf.map { meta, vcf -> vcf }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && scramble_mode == "yes" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(scramble_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf=merge_shortread_sv_callers.out
	}
	else if ( cnvnator_mode != "no" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf=merge_shortread_sv_callers.out
	}
	else if ( scramble_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(scramble_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf=merge_shortread_sv_callers.out
	}
	else {
		sv_vcf=Manta.out
	}

	sv_annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR_SV(sv_vcf,out_prefix,Phen2gene.out,sv_annovar_bed)
	SURVIVOR(ANNOVAR_SV.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)
	ngs_prio(out_prefix,RankVar.out,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,PhenoSV.out,ANNOVAR_SV.out,ANNOVAR.out.vcf_output,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}
