
process Manta {
	container ='beoungl/docker_test:manta'


        input:
        tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)


	output:
	path "${out_prefix}_manta.vcf"

	script:

	"""

	/manta/bin/configManta.py --bam=$bam --referenceFasta=$ref_fa --runDir ${out_prefix}_manta
	
	${out_prefix}_manta/runWorkflow.py -j 4 -g 64 -m local
	
	gunzip ${out_prefix}_manta/results/variants/diploidSV.vcf.gz

	mv ${out_prefix}_manta/results/variants/diploidSV.vcf ${out_prefix}_manta.vcf

	"""


}



