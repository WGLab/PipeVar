// Batch xTEA wrapper keyed by sample prefix for short-read mobile-element insertion calling.
process multi_xtea {
	container = 'beoungl/docker_test:xTea'

	input:
	tuple val(meta), path(bam), path(index)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(meta), path("${meta.id}_xtea.vcf"), emit: vcf
	tuple val(meta), path("${meta.id}_xtea_work"), optional: true, emit: work
	path "versions.yml", emit: versions

	script:
	def args = task.ext.args ?: ''
	def repLib = task.ext.rep_lib ?: '/opt/xtea/rep_lib_annotation'
	def gencode = task.ext.gencode_gff3 ?: '/opt/xtea/gencode.gff3'
	def xteaScripts = task.ext.xtea_scripts ?: '/opt/xTea/xtea'
	def xteaCmd = task.ext.xtea_cmd ?: 'xtea'
	def bamDotIndex = bam.name.endsWith('.cram') ? "${bam.name}.crai" : "${bam.name}.bai"
	def bamBaseIndex = bam.name.endsWith('.cram') ? bam.name.replaceFirst(/\.cram$/, '.crai') : bam.name.replaceFirst(/\.bam$/, '.bai')
	"""
	set -euo pipefail

	XTEA_WORK="${meta.id}_xtea_work"
	mkdir -p "\$XTEA_WORK"

	if [[ ! -d "${repLib}" ]]; then
	    echo "xTEA repeat library directory not found at ${repLib}" >&2
	    exit 1
	fi
	if [[ ! -s "${gencode}" ]]; then
	    echo "xTEA GENCODE GFF3 not found at ${gencode}" >&2
	    exit 1
	fi
	if [[ ! -d "${xteaScripts}" ]]; then
	    echo "xTEA script directory not found at ${xteaScripts}" >&2
	    exit 1
	fi

	if [[ "$bam" == *.cram ]]; then
	    [[ "$index" == "${bamDotIndex}" ]] || ln -sf "$index" "${bamDotIndex}"
	    [[ "$index" == "${bamBaseIndex}" ]] || ln -sf "$index" "${bamBaseIndex}"
	else
	    [[ "$index" == "${bamDotIndex}" ]] || ln -sf "$index" "${bamDotIndex}"
	    [[ "$index" == "${bamBaseIndex}" ]] || ln -sf "$index" "${bamBaseIndex}"
	fi

	printf '%s\\n' "${meta.id}" > sample_id.txt
	printf '%s\\t%s\\n' "${meta.id}" "$bam" > illumina_bam_list.txt

	set +u
	${xteaCmd} \\
	    -i sample_id.txt \\
	    -b illumina_bam_list.txt \\
	    -x null \\
	    -p "\$XTEA_WORK" \\
	    -o submit_jobs.sh \\
	    -l "${repLib}" \\
	    -r "$ref_fa" \\
	    -g "${gencode}" \\
	    --xtea "${xteaScripts}" \\
	    -f 5907 \\
	    -y 7 \\
	    --slurm \\
	    -t 0-72:00 \\
	    -q local \\
	    -n ${task.cpus} \\
	    -m 25 \\
	    $args
	set -u

	run_script=\$(find "\$XTEA_WORK" -type f -name 'run_xTea_pipeline.sh' | sort | head -n 1)
	if [[ -z "\$run_script" ]]; then
	    echo "xTEA did not generate run_xTea_pipeline.sh under \$XTEA_WORK" >&2
	    exit 1
	fi
	set +u
	bash "\$run_script"
	set -u

	xtea_vcf=\$(find "\$XTEA_WORK" -type f \\( -name '*.vcf' -o -name '*.gvcf' -o -name '*.vcf.gz' -o -name '*.gvcf.gz' \\) | sort | head -n 1)
	if [[ -n "\$xtea_vcf" ]]; then
	    if [[ "\$xtea_vcf" == *.gz ]]; then
	        gzip -cd "\$xtea_vcf" > ${meta.id}_xtea.vcf
	    else
	        cp "\$xtea_vcf" ${meta.id}_xtea.vcf
	    fi
	fi

	if [[ ! -s ${meta.id}_xtea.vcf ]]; then
	    {
	        echo '##fileformat=VCFv4.2'
	        echo "##source=xTEA"
	        echo "##reference=$ref_fa"
	        echo -e '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO'
	    } > ${meta.id}_xtea.vcf
	fi

	if ! grep -q '^#CHROM' ${meta.id}_xtea.vcf; then
	    echo "xTEA output is not a valid VCF: missing #CHROM header" >&2
	    exit 1
	fi

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  xtea: "unknown"
	END_VERSIONS
	"""

	stub:
	"""
	mkdir -p ${meta.id}_xtea_work
	{
	    echo '##fileformat=VCFv4.2'
	    echo "##source=xTEA-stub"
	    echo -e '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO'
	} > ${meta.id}_xtea.vcf
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  xtea: "stub"
	END_VERSIONS
	"""
}
