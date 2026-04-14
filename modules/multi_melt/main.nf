// Batch mobile-element insertion calling for short-read alignments using MELT Single mode.
process multi_melt {
	container = 'beoungl/docker_test:melt_2.2.2'

	input:
	tuple val(out_prefix), path(bam), path(index)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(out_prefix), path("${out_prefix}_melt.vcf")

	script:
	def args = task.ext.args ?: ''
	def meltResourceRoot = params.genome == 'grch38' ? '/opt/melt/resources/grch38' : '/opt/melt/resources/hg38'
	"""
	MELT_JAR=/opt/melt/MELT.jar
	MELT_GENE_BED=${meltResourceRoot}/genes.bed
	MELT_ZIP_LIST=${meltResourceRoot}/mei_zip_paths.txt

	if [[ ! -s "\$MELT_JAR" ]]; then
	    echo "MELT jar not found at \$MELT_JAR" >&2
	    exit 1
	fi
	if [[ ! -s "\$MELT_GENE_BED" ]]; then
	    echo "MELT gene BED not found at \$MELT_GENE_BED" >&2
	    exit 1
	fi
	if [[ ! -s "\$MELT_ZIP_LIST" ]]; then
	    echo "MELT zip list not found at \$MELT_ZIP_LIST" >&2
	    exit 1
	fi

	mapfile -t MELT_ZIPS < <(grep -v '^[[:space:]]*#' "\$MELT_ZIP_LIST" | awk 'NF')
	if [[ \${#MELT_ZIPS[@]} -eq 0 ]]; then
	    echo "No MELT transposon ZIPs configured in \$MELT_ZIP_LIST" >&2
	    exit 1
	fi

	declare -a FINAL_VCFS=()
	for zip_entry in "\${MELT_ZIPS[@]}"; do
	    if [[ "\$zip_entry" = /* ]]; then
	        zip_path="\$zip_entry"
	    else
	        zip_path="${meltResourceRoot}/\$zip_entry"
	    fi
	    if [[ ! -s "\$zip_path" ]]; then
	        echo "MELT transposon ZIP not found at \$zip_path" >&2
	        exit 1
	    fi

	    family=\$(basename "\$zip_path")
	    family=\${family%_MELT.zip}
	    family=\${family%.zip}
	    family_workdir="${out_prefix}_melt_runs/\${family}"
	    mkdir -p "\$family_workdir"

	    java -Xmx16G -jar "\$MELT_JAR" Single \\
	        -bamfile $bam \\
	        -h $ref_fa \\
	        -n "\$MELT_GENE_BED" \\
	        -t "\$zip_path" \\
	        -w "\$family_workdir" \\
	        $args

	    family_vcf=\$(find "\$family_workdir" -maxdepth 2 -name '*.final_comp.vcf' | head -n 1)
	    if [[ -z "\$family_vcf" ]]; then
	        echo "MELT did not emit a final_comp.vcf for \$family" >&2
	        exit 1
	    fi
	    FINAL_VCFS+=("\$family_vcf")
	done

	first_vcf="\${FINAL_VCFS[0]}"
	{
	    grep '^##' "\$first_vcf"
	    for vcf in "\${FINAL_VCFS[@]:1}"; do
	        grep '^##' "\$vcf" | grep -v '^##fileformat=' || true
	    done
	    grep '^#CHROM' "\$first_vcf" | tail -n 1
	    for vcf in "\${FINAL_VCFS[@]}"; do
	        grep -v '^#' "\$vcf"
	    done
	} > ${out_prefix}_melt.vcf
	"""
}
