// Batch xTEA wrapper keyed by sample prefix for short-read mobile-element insertion calling.
process multi_xtea {
	container = 'beoungl/docker_test:xtea_0.2'

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

	TASK_DIR="\$PWD"
	XTEA_WORK_NAME="${meta.id}_xtea_work"
	XTEA_NATIVE_WORK_NAME="${meta.id}_work"
	XTEA_WORK="\$TASK_DIR/\$XTEA_WORK_NAME"
	XTEA_NATIVE_WORK="\$TASK_DIR/\$XTEA_NATIVE_WORK_NAME"
	BAM_ABS="\$TASK_DIR/$bam"
	INDEX_ABS="\$TASK_DIR/$index"
	REF_ABS="\$TASK_DIR/$ref_fa"
	FA_INDEX_ABS="\$TASK_DIR/$fa_index"
	mkdir -p "\$XTEA_WORK"

	check_readable() {
	    local label="\$1"
	    local file_path="\$2"
	    if [[ ! -r "\$file_path" ]]; then
	        echo "xTEA input not readable: \${label} at \${file_path}" >&2
	        exit 1
	    fi
	}

	dump_xtea_context() {
	    echo "xTEA diagnostic context:" >&2
	    echo "  PWD=\$PWD" >&2
	    echo "  TASK_DIR=\$TASK_DIR" >&2
	    echo "  XTEA_WORK=\$XTEA_WORK" >&2
	    echo "  XTEA_NATIVE_WORK=\$XTEA_NATIVE_WORK" >&2
	    if [[ -e submit_jobs.sh ]]; then
	        echo "  submit_jobs.sh=present" >&2
	    else
	        echo "  submit_jobs.sh=absent" >&2
	    fi
	    echo "  TASK_DIR listing:" >&2
	    find "\$TASK_DIR" -maxdepth 2 -mindepth 1 -printf '    %p\\n' 2>/dev/null | sort | head -n 80 >&2 || true
	    echo "  XTEA_WORK listing:" >&2
	    find "\$XTEA_WORK" -maxdepth 4 -mindepth 1 -printf '    %p\\n' 2>/dev/null | sort | head -n 120 >&2 || true
	    echo "  XTEA_NATIVE_WORK listing:" >&2
	    find "\$XTEA_NATIVE_WORK" -maxdepth 4 -mindepth 1 -printf '    %p\\n' 2>/dev/null | sort | head -n 120 >&2 || true
	    if [[ -s xtea.generate.log ]]; then
	        echo "  xtea.generate.log tail:" >&2
	        tail -n 80 xtea.generate.log >&2
	    fi
	    if [[ -s xtea.run.log ]]; then
	        echo "  xtea.run.log tail:" >&2
	        tail -n 120 xtea.run.log >&2
	    fi
	    echo "  nested xTEA log tails:" >&2
	    find "\$XTEA_WORK" "\$XTEA_NATIVE_WORK" -type f \\( -name '*.log' -o -name '*.err' -o -name '*.out' -o -name '*.std_out' \\) 2>/dev/null \
	        | sort \
	        | head -n 40 \
	        | while IFS= read -r log_file; do
	            echo "    >>> \$log_file" >&2
	            tail -n 40 "\$log_file" >&2 || true
	          done
	}

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
	    [[ "$index" == "${bamDotIndex}" ]] || ln -sf "\$INDEX_ABS" "${bamDotIndex}"
	    [[ "$index" == "${bamBaseIndex}" ]] || ln -sf "\$INDEX_ABS" "${bamBaseIndex}"
	else
	    [[ "$index" == "${bamDotIndex}" ]] || ln -sf "\$INDEX_ABS" "${bamDotIndex}"
	    [[ "$index" == "${bamBaseIndex}" ]] || ln -sf "\$INDEX_ABS" "${bamBaseIndex}"
	fi

	check_readable "BAM/CRAM" "\$BAM_ABS"
	check_readable "BAM/CRAM index" "\$INDEX_ABS"
	check_readable "reference FASTA" "\$REF_ABS"
	check_readable "reference FASTA index" "\$FA_INDEX_ABS"
	export PYTHONFAULTHANDLER=1

	printf '%s\\n' "${meta.id}" > sample_id.txt
	printf '%s\\t%s\\n' "${meta.id}" "\$BAM_ABS" > illumina_bam_list.txt

	set +u
	if ! ${xteaCmd} \\
	    -i sample_id.txt \\
	    -b illumina_bam_list.txt \\
	    -x null \\
	    -p "\$XTEA_WORK" \\
	    -o submit_jobs.sh \\
	    -l "${repLib}" \\
	    -r "\$REF_ABS" \\
	    -g "${gencode}" \\
	    --xtea "${xteaScripts}" \\
	    -f 4883 \\
	    -y 7 \\
	    --slurm \\
	    -t 0-72:00 \\
	    -q local \\
	    -n ${task.cpus} \\
	    -m ${task.memory.toGiga()} \\
	    $args > xtea.generate.log 2>&1
	then
	    set -u
	    echo "xTEA command failed while generating run_xTEA_pipeline.sh" >&2
	    dump_xtea_context
	    exit 1
	fi
	set -u

	mapfile -t run_scripts < <(
	    {
	        if [[ -s submit_jobs.sh ]]; then
	            grep -o '[^[:space:]]*run_xTEA_pipeline[.]sh' submit_jobs.sh || true
	        fi
	        for search_root in "\$XTEA_WORK" "\$XTEA_NATIVE_WORK" "\$TASK_DIR"; do
	            [[ -e "\$search_root" ]] && find "\$search_root" -name 'run_xTEA_pipeline.sh' 2>/dev/null
	        done
	    } | sort -u
	)
	if (( \${#run_scripts[@]} == 0 )); then
	    echo "xTEA did not generate any run_xTEA_pipeline.sh under \$XTEA_WORK, \$XTEA_NATIVE_WORK, or \$TASK_DIR" >&2
	    dump_xtea_context
	    exit 1
	fi
	printf 'xTEA run scripts discovered:\\n' > xtea.run.log
	printf '  %s\\n' "\${run_scripts[@]}" >> xtea.run.log
	printf 'xTEA gVCF step bypassed; using -f 4883\\n' >> xtea.run.log

	set +u
	for run_script in "\${run_scripts[@]}"; do
	    run_dir="\$(dirname "\$run_script")"
	    run_name="\${run_script##*/}"
	    {
	        echo
	        echo "Running xTEA generated script: \$run_script"
	    } >> xtea.run.log
	    if ! ( cd "\$run_dir" && bash -Ee -o pipefail -x "./\$run_name" ) >> xtea.run.log 2>&1; then
	        set -u
	        echo "xTEA generated run script failed: \$run_script" >&2
	        dump_xtea_context
	        exit 1
	    fi
	done
	set -u

	xtea_candidates_to_vcf.py \\
	    --sample-id "${meta.id}" \\
	    --output ${meta.id}_xtea.vcf \\
	    --reference "\$REF_ABS" \\
	    --log xtea.run.log \\
	    "\$XTEA_WORK" "\$XTEA_NATIVE_WORK"

	if ! grep -q '^#CHROM' ${meta.id}_xtea.vcf; then
	    echo "xTEA output is not a valid VCF: missing #CHROM header" >&2
	    exit 1
	fi

	if [[ -d "\$XTEA_NATIVE_WORK_NAME" && -d "\$XTEA_WORK_NAME" ]] && ! find "\$XTEA_WORK_NAME" -mindepth 1 -print -quit | grep -q .; then
	    rmdir "\$XTEA_WORK_NAME"
	    ln -s "\$XTEA_NATIVE_WORK_NAME" "\$XTEA_WORK_NAME"
	elif [[ ! -e "\$XTEA_WORK_NAME" && -e "\$XTEA_NATIVE_WORK_NAME" ]]; then
	    ln -s "\$XTEA_NATIVE_WORK_NAME" "\$XTEA_WORK_NAME"
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
