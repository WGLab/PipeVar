// Merge and deduplicate short-read SV/MEI caller VCFs with Truvari collapse.
process truvari_shortread_sv_merge {
	container = 'beoungl/docker_test:truvari_0.1'

	input:
	path sv_vcfs
	tuple path(ref_fa), path(fa_index)
	val out_prefix

	output:
	path "${out_prefix}.shortread_sv.merged.vcf", emit: merged_vcf
	path "${out_prefix}.shortread_sv.truvari_collapsed.vcf", emit: collapsed_vcf
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

	mkdir -p truvari_inputs
	prepared_vcfs=()
	caller_names=(manta cnvnator xtea)
	manta_vcf=""

	for vcf in "\${vcf_files[@]}"; do
	    if [[ "\$(basename "\$vcf")" == *manta* ]]; then
	        manta_vcf="\$vcf"
	        break
	    fi
	done

	if [[ -z "\$manta_vcf" ]]; then
	    echo "Manta VCF is required to prepare short-read SV headers" >&2
	    exit 1
	fi

	if [[ "\$manta_vcf" == *.gz ]]; then
	    gzip -cd "\$manta_vcf"
	else
	    cat "\$manta_vcf"
	fi | awk '/^##/ { print; next } /^#CHROM/ { print; exit }' > truvari_inputs/manta.header.vcf

	for idx in "\${!vcf_files[@]}"; do
	    vcf="\${vcf_files[\$idx]}"
	    base_name="\$(basename "\$vcf")"
	    case "\$base_name" in
	        *xtea*) caller="xtea" ;;
	        *cnvnator*) caller="cnvnator" ;;
	        *manta*) caller="manta" ;;
	        *) caller="\${caller_names[\$idx]:-caller\$((idx + 1))}" ;;
	    esac
	    sample="${out_prefix}_\${caller}"
	    normalized="truvari_inputs/\${caller}.normalized.vcf"
	    prepared="truvari_inputs/\${caller}.prepared.vcf.gz"

	    if [[ "\$caller" == "cnvnator" ]]; then
	        cat truvari_inputs/manta.header.vcf
	        if [[ "\$vcf" == *.gz ]]; then
	            gzip -cd "\$vcf"
	        else
	            cat "\$vcf"
	        fi | awk '!/^#/'
	    elif [[ "\$vcf" == *.gz ]]; then
	        gzip -cd "\$vcf"
	    else
	        cat "\$vcf"
	    fi | awk -v caller="\$caller" -v sample="\$sample" '
	        BEGIN { OFS="\\t"; has_gt=0 }
	        /^##FORMAT=<ID=GT,/ { has_gt=1; print; next }
	        /^##/ { print; next }
	        /^#CHROM/ {
	            if (!has_gt) {
	                print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"
	            }
	            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"FORMAT",sample
	            next
	        }
	        NF >= 8 {
	            if (\$3 == "." || \$3 == "") {
	                \$3 = caller "_" NR
	            } else {
	                \$3 = caller "_" \$3
	            }
	            if (NF >= 10) {
	                print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10
	            } else {
	                print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"GT","0/1"
	            }
	        }
	    ' > "\$normalized"

	    bcftools sort -Oz -o "\$prepared" "\$normalized"
	    tabix -f -p vcf "\$prepared"
	    prepared_vcfs+=( "\$prepared" )
	done

	bcftools merge -m id "\${prepared_vcfs[@]}" -Oz -o ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz
	tabix -f -p vcf ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz

	truvari collapse \\
	    -i ${out_prefix}.shortread_sv.bcftools_merged.vcf.gz \\
	    -o ${out_prefix}.shortread_sv.truvari_merged.vcf \\
	    -c ${out_prefix}.shortread_sv.truvari_collapsed.vcf \\
	    -f $ref_fa \\
	    --intra \\
	    $args

	cp ${out_prefix}.shortread_sv.truvari_merged.vcf ${out_prefix}.shortread_sv.merged.vcf

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
