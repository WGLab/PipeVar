
def helpMessage = """
================================================================================
  P i p e V a r   P i p e l i n e
================================================================================
Usage:
  nextflow run main.nf --input_csv samples.csv [options]

--------------------------------------------------------------------------------
  MANDATORY INPUTS (Choose One Method)
--------------------------------------------------------------------------------
  --input_csv <path>       CSV manifest file containing sample information.
                           Columns required: sample, file_path, note_path
                           (Must use with --bam true OR --vcf true)

  --bam <path>             Path to a single BAM or CRAM file.
  
  --vcf <path>             Path to a single VCF file.

  --ref_fa <path>          Path to the reference FASTA file. 
                           (Required for all alignment/calling steps)

--------------------------------------------------------------------------------
  ANALYSIS MODE (Optional)
--------------------------------------------------------------------------------
  --mode <string>          Specific analysis to run.
                           Available: 'snp', 'sv'
                           Default: Runs BOTH SNP and SV analysis if omitted.
  
  --type <string>          Sequencing data type.
                           Available: 'ont' (default), 'pacbio', 'short'

--------------------------------------------------------------------------------
  FILTERING & TUNING
--------------------------------------------------------------------------------
  --genome <string>        Genome build (hg38, grch38). Default: hg38
  --light <yes|no>         'Light' mode (faster/lighter models). Default: no
  --gnomad <float>         Max gnomAD frequency. Default: 0.0001
  --rankscore <float>      Min RankScore. Default: 0.50
  --gq <int>               Min Genotype Quality. Default: 20
  --ad <int>               Min Allelic Depth. Default: 15
  --phen2gene_filter <int> Phen2Gene score threshold. Default: 500
  
================================================================================
"""


if (params.help) {
    println(helpMessage)
    exit(0)
}

// ------------------------------------------------------------------
// 0. PARAMETER NORMALIZATION (Fixing Capitalization/Spaces)
// ------------------------------------------------------------------

// We use the Elvis operator (?:) to handle nulls safely, defaulting to empty string if null
// We apply .trim() to remove spaces " ont " -> "ont"
// We apply .toLowerCase() to fix caps "ONT" -> "ont"

def clean_mode   = params.mode   ? params.mode.trim().toLowerCase()   : null
def clean_type   = params.type   ? params.type.trim().toLowerCase()   : 'ont' // Default from params
def clean_light  = params.light  ? params.light.trim().toLowerCase()  : 'no'
def clean_genome = params.genome ? params.genome.trim().toLowerCase() : 'hg38'
def clean_target = params.target ? params.target.trim().toLowerCase() : 'no'
def clean_note   = params.note   ? params.note.trim().toLowerCase()   : 'no'

// ------------------------------------------------------------------
// 1. INPUT VALIDATION (Catching Typos)
// ------------------------------------------------------------------

// Define allowed vocabularies
def valid_modes   = ['snp', 'sv']
def valid_types   = ['ont', 'pacbio', 'short']
def valid_genomes = ['hg38', 'grch38']
def valid_bools   = ['yes', 'no', 'true', 'false', null] // Allow null for safety

// CHECK 1: Validate Mode
if (clean_mode && !valid_modes.contains(clean_mode)) {
    error """
    ================================================================
    ERROR: Invalid Mode Specified
    ================================================================
    You provided: --mode "${params.mode}"
    
    We detected a typo. The valid options are:
      --mode snp
      --mode sv
      
    (Case insensitive, but spelling must match)
    ================================================================
    """
}

// CHECK 2: Validate Data Type
if (!valid_types.contains(clean_type)) {
    error """
    ================================================================
    ERROR: Invalid Data Type
    ================================================================
    You provided: --type "${params.type}"
    
    Valid options are:
      --type ont
      --type pacbio
      --type short
    ================================================================
    """
}

// CHECK 3: Validate Genome
if (!valid_genomes.contains(clean_genome)) {
    error """
    ================================================================
    ERROR: Unsupported Genome Build
    ================================================================
    You provided: --genome "${params.genome}"
    
    This pipeline currently supports:
      --genome hg38
      --genome grch38
    ================================================================
    """
}
// ------------------------------------------------------------------
// INPUT VALIDATION AND ERROR MESSAGES
// ------------------------------------------------------------------

// 1. Check for Missing Inputs
if (!params.input_csv && !params.bam && !params.vcf) {
    error """
    ERROR: No input data specified.
    You must provide one of the following:
      --input_csv <path>  : Path to CSV manifest
      --bam <path>        : Path to single BAM/CRAM file
      --vcf <path>        : Path to single VCF file
    """
}

// 2. Check for Ambiguous CSV Logic (Silent Failure Prevention)
// Your current logic requires params.bam or params.vcf to be truthy to process a CSV.
if (params.input_csv && !params.bam && !params.vcf) {
    error """
    ERROR: Input CSV specified, but data type is undefined.
    When using --input_csv, you must also specify the data type flag to trigger the correct logic:
      --bam true   (if CSV contains BAM/CRAM paths)
      --vcf true   (if CSV contains VCF paths)
    """
}

// 3. Check for Missing Reference (Mandatory for BAM/CRAM)
// Logic: If BAM is being processed (either single or via CSV), Ref is required.
def is_bam_mode = params.bam || (params.input_csv && params.bam)
if (is_bam_mode && !params.ref_fa) {
    error """
    ERROR: Reference FASTA missing.
    Processing BAM/CRAM files requires a reference genome.
    Please specify: --ref_fa <path/to/reference.fasta>
    """
}


// 5. Check for Ref Fasta Index (Optional but recommended)
// If ref is provided, ensure .fai exists to prevent downstream crashes
if (params.ref_fa) {
    def ref_fai = file("${params.ref_fa}.fai")
    if (!ref_fai.exists()) {
        error """
        ERROR: Reference index (.fai) not found.
        Expected at: ${ref_fai}
        Please index your reference: samtools faidx ${params.ref_fa}
        """
    }
}

include { SINGLE_ALIGNMENT_ALL_LONGPHASE } from './subworkflows/single_alignment_all_longphase' 
include { SINGLE_ALIGNMENT_ALL_NGS } from './subworkflows/single_alignment_all_ngs'
include { SINGLE_ALIGNMENT_ALL_LIGHT_LONGPHASE } from './subworkflows/single_alignment_all_light_longphase'
include { SINGLE_ALIGNMENT_ALL_NGS_LIGHT } from './subworkflows/single_alignment_all_ngs_light'
include { SINGLE_ALIGNMENT_LONG_SNP } from './subworkflows/single_alignment_long_snp'
include { SINGLE_ALIGNMENT_LONG_SV } from './subworkflows/single_alignment_long_sv'
include { SINGLE_ALIGNMENT_LONG_SNP_LIGHT } from './subworkflows/single_alignment_long_snp_light'
include { SINGLE_ALIGNMENT_NGS_SNP } from './subworkflows/single_alignment_ngs_snp'
include { SINGLE_ALIGNMENT_NGS_SV } from './subworkflows/single_alignment_ngs_sv'
include { SINGLE_ALIGNMENT_NGS_SNP_LIGHT } from './subworkflows/single_alignment_ngs_snp_light'
include { SINGLE_ALIGNMENT_VCF_SNP } from './subworkflows/single_alignment_vcf_snp'
include { SINGLE_ALIGNMENT_VCF_SV } from './subworkflows/single_alignment_vcf_sv'

include { INPUT_CSV_ALIGNMENT_ALL_LONGPHASE } from './subworkflows/input_csv_alignment_all_longphase'
include { INPUT_CSV_ALIGNMENT_ALL_LIGHT_LONGPHASE } from './subworkflows/input_csv_alignment_all_light_longphase'
include { INPUT_CSV_ALIGNMENT_ALL_NGS } from './subworkflows/input_csv_alignment_all_ngs'
include { INPUT_CSV_ALIGNMENT_ALL_NGS_LIGHT } from './subworkflows/input_csv_alignment_all_ngs_light'
include { INPUT_CSV_ALIGNMENT_LONG_SNP_LIGHT } from './subworkflows/input_csv_alignment_long_snp_light'
include { INPUT_CSV_ALIGNMENT_LONG_SNP } from './subworkflows/input_csv_alignment_long_snp'
include { INPUT_CSV_ALIGNMENT_LONG_SV } from './subworkflows/input_csv_alignment_long_sv'
include { INPUT_CSV_ALIGNMENT_NGS_SNP_LIGHT } from './subworkflows/input_csv_alignment_ngs_snp_light'
include { INPUT_CSV_ALIGNMENT_NGS_SV } from './subworkflows/input_csv_alignment_ngs_sv'
include { INPUT_CSV_NGS_SNP } from './subworkflows/input_csv_alignment_ngs_snp'
include { INPUT_CSV_ALIGNMENT_VCF_SNP } from './subworkflows/input_csv_alignment_vcf_snp'
include { INPUT_CSV_ALIGNMENT_VCF_SV } from './subworkflows/input_csv_alignment_vcf_sv'




workflow {
	if ( params.input_csv ) {
		if ( params.vcf ) {
		input_vcf = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def vcf_file = file(row.file_path, checkIfExists: true)
        def note_file = file(row.note_path, checkIfExists: true)
	def out_prefix = row.sample

        // Return a List/tuple in a specific order
        return tuple(
            row.sample,
            vcf_file,
            note_file,
        )
    			}
	input_bam=null
		}
		else if ( params.bam ) { 
input_bam = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def bam_file = file(row.file_path, checkIfExists: true)
        
        // FIX: Manipulate the full string path, avoiding 'parent' nullability issues
        // and handle both .bam.bai and .bai conventions automatically.
        def bai_path = bam_file.toString().replaceFirst(/\.bam$/, ".bai")
        def bai = file(bai_path)
        
        // If the standard .bai doesn't exist, try the .bam.bai convention
        if (!bai.exists()) {
             bai = file("${bam_file}.bai")
        }

        // CRAM override
        if (bam_file.name.endsWith('.cram')) {
            bai = file("${bam_file}.crai")
        }
        
        // Ensure index actually exists before returning tuple to avoid downstream staging errors
        if (!bai.exists()) {
            error "Index file not found for ${bam_file}. Looked for: ${bai_path} or ${bam_file}.bai"
        }

        def note_file = file(row.note_path, checkIfExists: true)
	def out_prefix = row.sample
	def bam_file_parent = bam_file.parent
	def note_file_parent = note_file.parent

        // Return the tuple
        return tuple(
            out_prefix,
            bam_file,
            bai,
            note_file,
        )
                        }
	input_vcf=null

		}
	}
        else if ( params.bam != null ) {
Channel
    .fromPath(params.bam ) 
    .map { file ->
        def meta = [ id: file.simpleName ]
        
        // Dynamic Logic: If it ends in cram, look for .crai, else .bai
        // Note: Check if your index naming is file.cram.crai or file.crai
        def index = file.name.endsWith('.cram') 
                    ? file.parent / "${file.name}.crai" 
                    : file.parent / "${file.name}.bai"
        
        return [ file, index ]
    }
    .set { bam }
                out_prefix=Channel.value(params.out_prefix)
        }
        else if ( params.vcf != null ) {
                vcf=Channel.value(params.vcf)
                out_prefix=Channel.value(params.out_prefix)
        }
	//Need to add dict file in case there is we are running haplotypecaller!!!!
	if (params.ref_fa != null && params.light == 'yes' && params.type == 'short' ) {
ref_fa = Channel
    .fromPath(params.ref_fa)
    .map { fa_file ->
        def fai_file = file("${fa_file}.fai")
        def dict_file = file("${fa_file.parent}/${fa_file.baseName}.dict")
        return [ fa_file, fai_file, dict_file ]
    }
    .first()
	}
        else if ( params.ref_fa != null ) {
ref_fa = Channel
    .fromPath(params.ref_fa)
    .map { fa_file ->
        def fai_file = file("${fa_file}.fai")
        return [ fa_file, fai_file ]
    }
    .first()
        }
	//Logic is bit confusing, but if there is note present, just use the note and store is_note to indiciate is_note for run phenotagger in later step.
	def is_note = "no"
	if ( params.note == 'yes' ) {
		is_note = "yes"
	}
	else if ( params.note != null && params.note != 'yes' ) {
		note=Channel.value(params.note)
		is_note = "yes"
	}
	//In case if there is hpo, then we assign note as HPO term directly, then skip phenotagger in later subworkflows. Use same variable name, but it's HPO in reality.
	if ( params.hpo != null && is_note == "no" ) { 
		note=Channel.value(params.hpo)
	}
	output_directory_check=file(params.output_directory)
        if ( !output_directory_check.exists() ) {
                output_directory_check.mkdirs()
        }
	type=Channel.value(params.type)
	gnomad=Channel.value(params.gnomad)
        rankscore_filter=Channel.value(params.rankscore)
	gq=Channel.value(params.gq)
	ad=Channel.value(params.ad)
	def target="null"
	if ( params.target == 'yes' ) {
		target = "yes"
	}
	if ( params.input_csv ) {
        if ( input_vcf != null ) {
            if ( params.mode == 'sv' ) {
                INPUT_CSV_ALIGNMENT_VCF_SV(input_vcf, ref_fa, is_note)
            }
            else if ( params.mode == 'snp' ) {
                INPUT_CSV_ALIGNMENT_VCF_SNP(input_vcf, ref_fa, rankscore_filter, gnomad, gq, ad, is_note)
            }
        }
        else if ( input_bam != null ) {
            if ( params.type == 'short' ) {
                if ( params.light == 'yes' ) {
                    if ( params.mode == 'snp' ) {
                        INPUT_CSV_ALIGNMENT_NGS_SNP_LIGHT(input_bam, ref_fa, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                    else {
                        INPUT_CSV_ALIGNMENT_ALL_NGS_LIGHT(input_bam, ref_fa, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                }
                else { // Not light
                    if ( params.mode == 'sv' ) {
                        INPUT_CSV_ALIGNMENT_NGS_SV(input_bam, ref_fa,  is_note)
                    }
                    else if ( params.mode == 'snp' ) {
                        INPUT_CSV_NGS_SNP(input_bam, ref_fa,  rankscore_filter, gnomad, gq, ad, is_note, target)
			}
			else {
			INPUT_CSV_ALIGNMENT_ALL_NGS(input_bam, ref_fa,  rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                }
            }
            else { // Long reads
                if ( params.light == 'yes' ) {
                    if ( params.mode == 'snp' ) {
                        INPUT_CSV_ALIGNMENT_LONG_SNP_LIGHT(input_bam, ref_fa,  rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                    else {
                        INPUT_CSV_ALIGNMENT_ALL_LIGHT_LONGPHASE(input_bam, ref_fa,  rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                } // You were missing this closing brace
                else { // Not light
                    if ( params.mode == 'sv' ) {
                        INPUT_CSV_ALIGNMENT_LONG_SV(input_bam, ref_fa,  is_note)
                    }
                    else if ( params.mode == 'snp' ) { // Fixed missing quote here
                        INPUT_CSV_ALIGNMENT_LONG_SNP(input_bam, ref_fa,  rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                    else {
                        INPUT_CSV_ALIGNMENT_ALL_LONGPHASE(input_bam, ref_fa,  rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                }
            }
        }
    }
    else { // Single File Mode
        if ( params.vcf ) {
            if ( params.mode == 'sv' ) {
                SINGLE_ALIGNMENT_VCF_SV(vcf, out_prefix, ref_fa,  note, is_note)
            }
            else if ( params.mode == 'snp' ) {
                SINGLE_ALIGNMENT_VCF_SNP(vcf, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note)
            }
        }
        else if ( params.bam != null ) {
            if ( params.type == 'short' ) {
                if ( params.light == 'yes' ) {
                    if ( params.mode == 'snp' ) {
                        SINGLE_ALIGNMENT_NGS_SNP_LIGHT(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                    else {
                        SINGLE_ALIGNMENT_ALL_NGS_LIGHT(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                }
                else {
                    if ( params.mode == 'sv' ) {
                        SINGLE_ALIGNMENT_NGS_SV(bam, out_prefix, ref_fa,  note, is_note)
                    }
                    else if ( params.mode == 'snp' ) {
                        SINGLE_ALIGNMENT_NGS_SNP(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
		else {
		SINGLE_ALIGNMENT_ALL_NGS(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
		}
                }
            }
            else { // Long reads
                if ( params.light == 'yes' ) {
                    if ( params.mode == 'snp' ) {
                        // Fixed is_not typo below
                        SINGLE_ALIGNMENT_LONG_SNP_LIGHT(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                    else {
                        SINGLE_ALIGNMENT_ALL_LIGHT_LONGPHASE(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                } // Missing brace fixed here
                else {
                    if ( params.mode == 'sv' ) {
                        SINGLE_ALIGNMENT_LONG_SV(bam, out_prefix, ref_fa,  note, is_note)
                    }
                    else if ( params.mode == 'snp' ) { // Fixed missing quote here
                        // Fixed is_not typo below
                        SINGLE_ALIGNMENT_LONG_SNP(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                    else {
                        SINGLE_ALIGNMENT_ALL_LONGPHASE(bam, out_prefix, ref_fa,  note, rankscore_filter, gnomad, gq, ad, is_note, target)
                    }
                }
            }
        }
    }
}



workflow.onComplete {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
        println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

