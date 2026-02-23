
// Run RankVar to prioritize small variants with phenotype-aware evidence.
process RankVar {
        container ='beoungl/docker_test:rankvar'
	
	input:
	path vcf 
	path phen2gene
	path hpo
	val out_prefix
	val gnomad
	val gq
	val ad
	val rankvar_filter

	output:
	path "${out_prefix}.rank_var.tsv"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate rankvar

	# Prefilter ANNOVAR table before RankVar to reduce input size:
	# exonic/splicing only + non-reference GT + GQ threshold from Otherinfo12/Otherinfo13.
	awk -F'\t' -v OFS='\t' -v min_gq="$gq" '
	function trim(x) { gsub(/^[ \t]+|[ \t]+\$/, "", x); return x }
	function pass_gt_gq(fmt, sample,   fmt_a,smp_a,nf,ns,i,gt_i,gq_i,gt,gqv) {
	    gt_i=0; gq_i=0; gt=""; gqv=""
	    nf=split(fmt, fmt_a, ":"); ns=split(sample, smp_a, ":")
	    for (i=1; i<=nf; i++) {
	        if (fmt_a[i] == "GT") gt_i=i
	        else if (fmt_a[i] == "GQ") gq_i=i
	    }
	    if (gt_i>0 && gt_i<=ns) gt=trim(smp_a[gt_i])
	    if (gq_i>0 && gq_i<=ns) gqv=trim(smp_a[gq_i])
	    if (gt=="" || gt=="./." || gt==".|." || gt=="0/0" || gt=="0|0") return 0
	    if (gqv!="" && gqv!="." && (gqv+0) < min_gq) return 0
	    return 1
	}
	NR==1 {
	    for (i=1; i<=NF; i++) h[\$i]=i
	    fmt_idx = (("Otherinfo12" in h) ? h["Otherinfo12"] : 0)
	    smp_idx = (("Otherinfo13" in h) ? h["Otherinfo13"] : 0)
	    print
	    next
	}
	(\$6=="exonic" || \$6=="splicing") {
	    pass_qc = 1
	    if (fmt_idx>0 && smp_idx>0) pass_qc = pass_gt_gq(\$(fmt_idx), \$(smp_idx))
	    if (pass_qc) print
	}' $vcf > ${out_prefix}.rankvar_temp.txt

	python /opt/RankVar/RankVar.py --annovar ${out_prefix}.rankvar_temp.txt --output ${out_prefix}_rankvar --hpo_ids $hpo --phen2gene $phen2gene --gq $gq --ad $ad --gnomad $gnomad
	
	#mv ${out_prefix}_rankvar/rank_var.tsv ${out_prefix}.rank_var.tsv
	# This pathogenicity score cutoff is configurable via workflow param.
	awk -F'\t' -v cutoff="$rankvar_filter" 'NR==1 || \$12 > cutoff' ${out_prefix}_rankvar/rank_var.tsv > ${out_prefix}.rank_var.tsv

	rm ${out_prefix}.rankvar_temp.txt
	"""
}
