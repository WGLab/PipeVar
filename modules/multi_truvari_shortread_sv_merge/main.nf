// Batch merge and deduplicate short-read SV/MEI caller VCFs with Truvari collapse.
process multi_truvari_shortread_sv_merge {
	container = 'beoungl/docker_test:truvari_0.2'

	input:
	tuple val(out_prefix), path(sv_vcfs)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(out_prefix), path("${out_prefix}.shortread_sv.merged.vcf"), emit: merged_vcf
	tuple val(out_prefix), path("${out_prefix}.shortread_sv.truvari_collapsed.vcf"), emit: collapsed_vcf
	path "versions.yml", emit: versions

	script:
	def args = task.ext.args ?: ''
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari

	vcf_files=( $sv_vcfs )
	if [[ \${#vcf_files[@]} -eq 0 ]]; then
	    echo "No short-read SV VCFs provided for Truvari merge" >&2
	    exit 1
	fi

	prepare_shortread_sv_merge.py \\
	    --out-prefix ${out_prefix} \\
	    --reference "$ref_fa" \\
	    --work-dir truvari_inputs \\
	    --truvari-args "$args" \\
	    "\${vcf_files[@]}"

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  truvari: "\$(truvari version 2>/dev/null || echo unknown)"
	  bcftools: "\$(bcftools --version | head -n 1 | awk '{print \$2}')"
	END_VERSIONS
	"""

	stub:
	"""
	{
	    echo '##fileformat=VCFv4.2'
	    echo '##source=truvari-collapse-stub'
	    echo -e '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO'
	} > ${out_prefix}.shortread_sv.merged.vcf
	cp ${out_prefix}.shortread_sv.merged.vcf ${out_prefix}.shortread_sv.truvari_collapsed.vcf
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  truvari: "stub"
	  bcftools: "stub"
	END_VERSIONS
	"""
}
