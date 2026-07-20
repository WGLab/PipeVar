// Batch merge and deduplicate short-read SV/MEI caller VCFs with Truvari collapse.
process multi_truvari_shortread_sv_merge {
	container = 'beoungl/docker_test:truvari_0.5'

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

	/usr/local/bin/prepare_shortread_sv_merge.py --version >&2
	/usr/local/bin/prepare_shortread_sv_merge.py \\
	    --out-prefix ${out_prefix} \\
	    --reference "$ref_fa" \\
	    --work-dir truvari_inputs \\
	    --prepare-only \\
	    "\${vcf_files[@]}"

	prepared_vcfs=()
	for normalized_vcf in truvari_inputs/*.normalized.vcf; do
	    if [[ ! -s "\$normalized_vcf" ]]; then
	        echo "Prepared normalized VCF is missing or empty: \$normalized_vcf" >&2
	        exit 1
	    fi
	    prepared_vcf="\${normalized_vcf%.normalized.vcf}.prepared.vcf.gz"
	    bcftools sort -Oz -o "\$prepared_vcf" "\$normalized_vcf"
	    tabix -f -p vcf "\$prepared_vcf"
	    prepared_vcfs+=( "\$prepared_vcf" )
	done

	if [[ \${#prepared_vcfs[@]} -eq 0 ]]; then
	    echo "No normalized short-read SV VCFs were prepared for Truvari merge" >&2
	    exit 1
	fi

	bcftools merge -m id "\${prepared_vcfs[@]}" -Oz -o ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz
	tabix -f -p vcf ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz

	bcftools view -e 'INFO/SVTYPE="BND"' ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz -Oz -o ${out_prefix}.shortread_sv.collapse_input.vcf.gz
	tabix -f -p vcf ${out_prefix}.shortread_sv.collapse_input.vcf.gz
	bcftools view -i 'INFO/SVTYPE="BND"' ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz -Oz -o ${out_prefix}.shortread_sv.bnd_passthrough.vcf.gz
	tabix -f -p vcf ${out_prefix}.shortread_sv.bnd_passthrough.vcf.gz

	if bcftools view -H ${out_prefix}.shortread_sv.collapse_input.vcf.gz | grep -q .; then
	    truvari collapse \\
	        -i ${out_prefix}.shortread_sv.collapse_input.vcf.gz \\
	        -o ${out_prefix}.shortread_sv.truvari_merged.vcf \\
	        -c ${out_prefix}.shortread_sv.truvari_collapsed.vcf \\
	        -f "$ref_fa" \\
	        --intra $args
	else
	    bcftools view ${out_prefix}.shortread_sv.collapse_input.vcf.gz > ${out_prefix}.shortread_sv.truvari_merged.vcf
	    cp ${out_prefix}.shortread_sv.truvari_merged.vcf ${out_prefix}.shortread_sv.truvari_collapsed.vcf
	fi

	validate_vcf_width() {
	    local vcf="\$1"
	    awk -F '\\t' '
	        /^#CHROM/ {
	            expected = NF
	            next
	        }
	        /^#/ { next }
	        {
	            if (!expected) {
	                printf "VCF width validation failed: missing #CHROM header\\n" > "/dev/stderr"
	                exit 1
	            }
	            if (NF != expected) {
	                printf "VCF width validation failed at line %d: expected %d fields, observed %d\\n", NR, expected, NF > "/dev/stderr"
	                exit 1
	            }
	        }
	        END {
	            if (!expected) {
	                exit 1
	            }
	        }
	    ' "\$vcf"
	}

	recombine_bnd_with_target_schema() {
	    local target_vcf="\$1"
	    local output_vcf="\$2"
	    local label="\$3"
	    local sample_file="\${label}.samples.txt"
	    local target_bgz="\${label}.target.vcf.gz"
	    local projected_bnd_unsorted="\${label}.bnd.projected.unsorted.vcf.gz"
	    local projected_bnd="\${label}.bnd.projected.vcf.gz"
	    local combined_unsorted="\${label}.combined.unsorted.vcf"

	    bcftools query -l "\$target_vcf" > "\$sample_file"
	    if [[ -s "\$sample_file" ]]; then
	        local sample_csv
	        sample_csv=\$(paste -sd, "\$sample_file")
	        bcftools view -s "\$sample_csv" -Oz -o "\$projected_bnd_unsorted" ${out_prefix}.shortread_sv.bnd_passthrough.vcf.gz
	    else
	        bcftools view -G -Oz -o "\$projected_bnd_unsorted" ${out_prefix}.shortread_sv.bnd_passthrough.vcf.gz
	    fi
	    bcftools sort -Oz -o "\$projected_bnd" "\$projected_bnd_unsorted"
	    tabix -f -p vcf "\$projected_bnd"

	    bcftools sort -Oz -o "\$target_bgz" "\$target_vcf"
	    tabix -f -p vcf "\$target_bgz"
	    bcftools concat -a -Ov "\$target_bgz" "\$projected_bnd" > "\$combined_unsorted"
	    bcftools sort -Ov -o "\$output_vcf" "\$combined_unsorted"
	    validate_vcf_width "\$output_vcf"
	}

	sort_plain_vcf() {
	    local input_vcf="\$1"
	    local output_vcf="\$2"
	    local sorted_vcf="\${output_vcf%.vcf}.sorted.tmp.vcf"

	    bcftools sort -Ov -o "\$sorted_vcf" "\$input_vcf"
	    mv "\$sorted_vcf" "\$output_vcf"
	}

	if bcftools view -H ${out_prefix}.shortread_sv.bnd_passthrough.vcf.gz | grep -q .; then
	    echo "Passing BND records through without Truvari collapse using each target VCF sample schema" >&2
	    recombine_bnd_with_target_schema \\
	        ${out_prefix}.shortread_sv.truvari_merged.vcf \\
	        ${out_prefix}.shortread_sv.merged.vcf \\
	        retained
	    recombine_bnd_with_target_schema \\
	        ${out_prefix}.shortread_sv.truvari_collapsed.vcf \\
	        ${out_prefix}.shortread_sv.truvari_collapsed.with_bnd.vcf \\
	        collapsed
	    mv ${out_prefix}.shortread_sv.truvari_collapsed.with_bnd.vcf ${out_prefix}.shortread_sv.truvari_collapsed.vcf
	else
	    sort_plain_vcf \\\\
	        ${out_prefix}.shortread_sv.truvari_merged.vcf \\\\
	        ${out_prefix}.shortread_sv.merged.vcf
	    sort_plain_vcf \\\\
	        ${out_prefix}.shortread_sv.truvari_collapsed.vcf \\\\
	        ${out_prefix}.shortread_sv.truvari_collapsed.vcf
	fi

	validate_vcf_width ${out_prefix}.shortread_sv.merged.vcf
	validate_vcf_width ${out_prefix}.shortread_sv.truvari_collapsed.vcf

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
