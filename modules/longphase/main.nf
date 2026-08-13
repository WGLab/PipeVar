
// Phase SNV/SV evidence and aggregate ranked evidence into final prioritized VCF.
process longphase {
	container ='beoungl/docker_test:longphase_0.2.35'

	input:
	tuple path(bam_path), path(index)
	path annovar_vcf
	path sv_phase_vcf
	path sv_annotation_vcf
	path sv_pathogenic
	path snv_rankscore
	path snv_pathogenic
	path snv_rankvar
	path hpo_path
	val out_prefix
	tuple path(ref_fa), path(fa_index)
	val(inheritance_mode)
	val(include_clinvar_report)
	val(allow_unphased_comphet)

	output:
	path "${out_prefix}.prio.vcf"
	path "${out_prefix}.prio_gene.vcf"

	
	script:
	def platform_args = task.ext.args != null ? task.ext.args : '--ont'
	def min_score = params.phenosv_score ?: '0.50'
	def gene_filter = params.gene ?: ''
	def sv_only = params.prioritize_sv_only ?: 'no'


		"""
		bcftools reheader --fai $fa_index --output ${out_prefix}.longphase_input.vcf $annovar_vcf

		awk '
		BEGIN {
		    while ((getline line < ARGV[1]) > 0) {
		        split(line, columns, "\t")
		        reference_length[columns[1]] = columns[2]
		    }
		    close(ARGV[1])
		    ARGV[1] = ""
		}
		/^##contig=</ {
		    contig_headers++
		    if (\$0 !~ /,length=[1-9][0-9]*([,>])/) {
		        print "ERROR: LongPhase SNV VCF contig header lacks a positive numeric length: " \$0 > "/dev/stderr"
		        errors = 1
		    }
		    next
		}
		/^#/ { next }
		{
		    if (!(\$1 in reference_length)) {
		        print "ERROR: LongPhase SNV VCF contains contig absent from reference FAI at " \$1 ":" \$2 > "/dev/stderr"
		        errors = 1
		    }
		    if (\$2 !~ /^[1-9][0-9]*\$/) {
		        print "ERROR: LongPhase SNV VCF contains an invalid POS at " \$1 ":" \$2 > "/dev/stderr"
		        errors = 1
		    }
		    else if ((\$1 in reference_length) && \$2 > reference_length[\$1]) {
		        print "ERROR: LongPhase SNV VCF POS exceeds reference contig length at " \$1 ":" \$2 > "/dev/stderr"
		        errors = 1
		    }
		}
		END {
		    if (contig_headers == 0) {
		        print "ERROR: LongPhase SNV VCF has no contig headers after FAI reheadering" > "/dev/stderr"
		        errors = 1
		    }
		    exit errors
		}' $fa_index ${out_prefix}.longphase_input.vcf

		/longphase_linux-x64 phase -s ${out_prefix}.longphase_input.vcf --sv-file=$sv_phase_vcf -t ${task.cpus} -o ${out_prefix}_phased $platform_args -b $bam_path -r $ref_fa

		/longphase_linux-x64 haplotag -r $ref_fa -s ${out_prefix}_phased.vcf --sv-file ${out_prefix}_phased_SV.vcf -b $bam_path -t ${task.cpus} -o ${out_prefix}_haplotag

	        bash /phenosv_vcf_and_tsv.sh $sv_pathogenic ${out_prefix}_phased_SV.vcf $sv_annotation_vcf $out_prefix
		mv ${out_prefix}.phenosv.vcf ${out_prefix}.phenosv.unfiltered.vcf
		python3 /filter_phenosv_vcf.py ${out_prefix}.phenosv.unfiltered.vcf ${out_prefix}.phenosv.vcf --min-score $min_score --genes "$gene_filter"

	if [[ "$sv_only" == "yes" ]]; then
	    python3 /assign_dom_or_rec_sv_only.py ${out_prefix}.phenosv.vcf $hpo_path $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"
	else
	    bash /clinvar_vcf_and_txt.sh $snv_pathogenic ${out_prefix}_phased.vcf $out_prefix
	    bash /rankscore_vcf_and_txt.sh $snv_rankscore ${out_prefix}_phased.vcf $out_prefix
	    bash /rankvar_vcf_and_tsv.sh $snv_rankvar ${out_prefix}_phased.vcf $out_prefix
	    python3 /assign_dom_or_rec.py ${out_prefix}.clinvar.vcf ${out_prefix}.phenosv.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf $hpo_path $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"
	fi

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio_gene.vcf gene --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter"
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf variant --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter"




	"""

}
