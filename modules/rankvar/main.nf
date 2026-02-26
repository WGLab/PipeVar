
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

	# Replace ANNOVAR-normalized coordinates with VCF-origin coordinates (Otherinfo4/5/7/8)
	# to preserve matching with phased VCF records in splicing variants.
	awk -F'\t' -v OFS='\t' '
	function fail(msg) { print "ERROR: " msg > "/dev/stderr"; exit 2 }
	FNR==NR {
	    if (FNR == 1) {
	        for (i=1; i<=NF; i++) h[\$i]=i
	        chr_i = h["Chr"]; start_i = h["Start"]; ref_i = h["Ref"]; alt_i = h["Alt"]
	        o4_i = h["Otherinfo4"]; o5_i = h["Otherinfo5"]; o7_i = h["Otherinfo7"]; o8_i = h["Otherinfo8"]
	        if (!chr_i || !start_i || !ref_i || !alt_i) fail("ANNOVAR TXT missing Chr/Start/Ref/Alt columns for RankVar remap")
	        if (!o4_i || !o5_i || !o7_i || !o8_i) fail("ANNOVAR TXT missing Otherinfo4/5/7/8 columns for RankVar remap")
	        next
	    }
	    key = \$(chr_i) SUBSEP \$(start_i) SUBSEP \$(ref_i) SUBSEP \$(alt_i)
	    val = \$(o4_i) SUBSEP \$(o5_i) SUBSEP \$(o7_i) SUBSEP \$(o8_i)
	    if ((key in map) && map[key] != val) fail("Conflicting ANNOVAR->VCF coordinate mapping for key: " key)
	    map[key] = val
	    next
	}
	FNR==1 {
	    delete h2
	    for (i=1; i<=NF; i++) h2[\$i]=i
	    r_chr = h2["Chr"]; r_start = h2["Start"]; r_ref = h2["Ref"]; r_alt = h2["Alt"]
	    if (!r_chr || !r_start || !r_ref || !r_alt) fail("rank_var.tsv missing required header(s): Chr, Start, Ref, Alt")
	    print
	    next
	}
	{
	    key = \$(r_chr) SUBSEP \$(r_start) SUBSEP \$(r_ref) SUBSEP \$(r_alt)
	    if (key in map) {
	        n = split(map[key], a, SUBSEP)
	        \$(r_chr)=a[1]; \$(r_start)=a[2]; \$(r_ref)=a[3]; \$(r_alt)=a[4]
	    }
	    print
	}
	' $vcf ${out_prefix}.rank_var.tsv > ${out_prefix}.rank_var.tsv.tmp && mv ${out_prefix}.rank_var.tsv.tmp ${out_prefix}.rank_var.tsv

	rm ${out_prefix}.rankvar_temp.txt
	"""
}
