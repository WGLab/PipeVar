
// Batch structural-variant calling from short-read BAM inputs with Manta.
process multi_manta {
	container ='beoungl/docker_test:manta'


        input:
        tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)


	output:
	tuple val(out_prefix), path("${out_prefix}_manta.vcf") 

	script:

	"""

	/manta/bin/configManta.py --bam=$bam --referenceFasta=$ref_fa --runDir ${out_prefix}_manta
	
	${out_prefix}_manta/runWorkflow.py -j ${task.cpus} -g ${task.memory.toGiga()} -m local
	
	gunzip ${out_prefix}_manta/results/variants/diploidSV.vcf.gz

	mv ${out_prefix}_manta/results/variants/diploidSV.vcf ${out_prefix}_manta.vcf

	"""


}


