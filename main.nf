def pipelineVersion = params.pipeline_version ?: "0.3.0"

def helpMessage = """
================================================================================
  P i p e V a r   P i p e l i n e
  Version: ${pipelineVersion}
================================================================================
QUICK USAGE
  Single sample mode:
    nextflow run main.nf [--bam sample.bam | --vcf sample.vcf] [options]
    nextflow run main.nf --annotated_snv yes --annovar_txt sample.hg38_multianno.txt --vcf sample.hg38_multianno.vcf [options]

  CSV batch mode:
    nextflow run main.nf --input_csv samples.csv --bam true [options]
    nextflow run main.nf --input_csv samples.csv --vcf true [options]

--------------------------------------------------------------------------------
EXECUTION PROFILES (nextflow.config)
--------------------------------------------------------------------------------
  -profile standard
      Default. SLURM executor + Singularity backend.

  -profile slurm_singularity
      Explicit SLURM + Singularity.

  -profile local_singularity
      Local executor + Singularity.

  -profile local_docker
      Local executor + Docker.

  Container mount paths (used by Singularity/Docker profiles):
    --annovar_host_path <DIR>   Host path mounted to /annovar
    --phenosv_host_path <DIR>   Host path mounted to /PhenoSV/train_data

--------------------------------------------------------------------------------
INPUT MODES
--------------------------------------------------------------------------------
  1) Single sample BAM/CRAM mode
     Required:
       --bam <FILE>
       --ref_fa <FILE>
       For short-read runs (`--type short`):
         --expansionhunter_variant_catalog <FILE>
       One phenotype source:
         --note <FILE>   (clinical note; pipeline runs PhenoTagger)
         --hpo  <FILE>   (HPO IDs; pipeline skips PhenoTagger)

  2) Single sample VCF mode
     Required:
       --vcf <FILE>
       --ref_fa <FILE>
       --mode <snp|sv>
       One phenotype source:
         --note <FILE> or --hpo <FILE>

  3) Single sample annotated-SNV mode
     Required:
       --annotated_snv yes
       --annovar_txt <FILE>
       --vcf <FILE>   (matching ANNOVAR multianno VCF)
       One phenotype source:
         --note <FILE> or --hpo <FILE>
     Notes:
       Without --bam: SNP prioritization only, --ref_fa not required
       With --bam/--type short: also runs short-read SV/CNV/STR and optional mito analysis
       With --bam/--type short:
         --expansionhunter_variant_catalog <FILE>
       --mode sv is not supported
       --target yes is not supported

  4) CSV batch BAM mode
     Required:
       --input_csv <FILE> --bam true --ref_fa <FILE>
       For short-read runs (`--type short`):
         --expansionhunter_variant_catalog <FILE>
     CSV columns:
       sample,file_path,note_path
       Optional age field for CSV prioritization flows:
         sample,file_path,note_path,age_of_onset
         sample,file_path,note_path,age
       Notes:
         age_of_onset is preferred if both columns exist
         age is interpreted per CSV row (per sample), not globally
         empty age is allowed (treated as not provided)
         non-empty age must be xd/xm/xy or integer years
         examples: 10d, 9m, 7y, 7 (normalized to 7y)
     Note handling:
       default            -> note_path treated as clinical note (PhenoTagger ON)
       --note no          -> note_path treated as HPO file (PhenoTagger OFF)

  5) CSV batch unified manifest mode
     Required:
       --input_csv <FILE>
     Required common columns:
       sample,input_kind,phenotype_path,phenotype_format
     Conditional columns:
       snv_txt_path,snv_vcf_path,vcf_path,alignment_path,alignment_index_path
     input_kind values:
       annotated_snv,vcf_snv,vcf_sv,bam_ngs,cram_ngs
     phenotype_format values:
       clinical_note,hpo
     Notes:
       annotated_snv rows may also include alignment_path/alignment_index_path to run SV/CNV/STR
       and optional mito analysis alongside imported SNP evidence
       annotated_snv rows with alignment_path and --type short also require:
         --expansionhunter_variant_catalog <FILE>

  6) CSV batch VCF mode
     Required:
       --input_csv <FILE> --vcf true --ref_fa <FILE> --mode <snp|sv>
     CSV columns:
       sample,file_path,note_path
       Optional age field for CSV prioritization flows:
         sample,file_path,note_path,age_of_onset
         sample,file_path,note_path,age
       Notes:
         age_of_onset is preferred if both columns exist
         age is interpreted per CSV row (per sample), not globally
         empty age is allowed (treated as not provided)
         non-empty age must be xd/xm/xy or integer years
         examples: 10d, 9m, 7y, 7 (normalized to 7y)

--------------------------------------------------------------------------------
CORE OPTIONS
--------------------------------------------------------------------------------
  --mode <snp|sv>          Run only SNP or SV branch. Omit to run both where applicable.
  --annotated_snv <yes|no> Use pre-annotated ANNOVAR TXT + VCF SNP input. Default: no
  --annovar_txt <FILE>     Single-sample ANNOVAR multianno TXT paired with --vcf
  --type <ont|pacbio|short>
                           Sequencing type for BAM/CRAM flows. Default: ont
  --light <yes|no>         Use lightweight SNP/SV models where supported. Default: no
  --expansionhunter_variant_catalog <FILE>
                           ExpansionHunter variant catalog for short-read BAM/CRAM workflows
  --genome <hg38|grch38>   Reserved for non-ExpansionHunter uses. Default: hg38

--------------------------------------------------------------------------------
FILTERING OPTIONS
--------------------------------------------------------------------------------
  --gnomad <FLOAT>         Max gnomAD AF filter for SNP prioritization. Default: 0.0001
  --inheritance_mode <ml|omim|gnomad>
                           Dominant/recessive assignment mode. Default: ml (ML-first, OMIM fallback)
                           Note: gnomad uses gnomAD constraint (LOEUF) fallback lists.
  --include_clinvar_report <yes|no>
                           Include ClinVar-only calls in final prioritized report outputs. Default: yes
  --allow_unphased_comphet <yes|no>
                           Treat unphased 0/1 or 1/0 AR het pairs as compound het in final prioritization. Default: no
  --rankscore <FLOAT>      Minimum RankScore cutoff. Default: 0.50
  --rankscore_softwares <CSV> Comma-separated RankScore software list. Default: all
  --rankvar <FLOAT>        Minimum RankVar score cutoff. Default: 0.05
  --gq <INT>               Minimum genotype quality. Default: 20
  --ad <INT>               Minimum allele depth. Default: 15
  --phen2gene_filter <INT> Number of top Phen2Gene genes used for targeted mode. Default: 500
  --target <yes|no>        If yes, run targeted calling in phenotype-derived gene regions.
  --cnvnator <yes|no>      Add CNVnator to short-read SV/all-NGS calling. Default: yes
  --cnvnator_bin_size <INT>
                           CNVnator read-depth bin size. Default: 100
  --scramble <yes|no>      Add SCRAMBLE mobile-element calling to short-read SV/all-NGS BAM/CRAM paths. Default: no
  --mito <yes|no>          Add mitochondrial Mutect2 + mtDNA annotation branch for short-read BAM/CRAM. Default: no
                           Requires BWA-indexed reference sidecars alongside --ref_fa.
  --mito_contig <STRING>   Preferred mitochondrial contig alias. Default: chrM
  --mito_min_vaf <FLOAT>   Minimum heteroplasmy/allele fraction retained in mito prioritization. Default: 0.01
  --mito_min_depth <INT>   Minimum depth retained in mito prioritization. Default: 50
  --mito_min_alt_reads <INT>
                           Minimum alternate read count retained in mito prioritization. Default: 5

--------------------------------------------------------------------------------
OUTPUT / GENERAL
--------------------------------------------------------------------------------
  --out_prefix <STRING>    Output prefix. Default: PipeVar
  --output_directory <DIR> Output directory. Default: launch directory
  --help                   Print this help text and exit

--------------------------------------------------------------------------------
EXAMPLES
--------------------------------------------------------------------------------
  1) SLURM + Singularity (default), single long-read BAM, full path:
     nextflow run main.nf \\
       -profile standard \\
       --bam /data/sample.bam \\
       --ref_fa /refs/hg38.fa \\
       --note /data/note.txt \\
       --out_prefix patient1 \\
       --type ont

  2) Local + Docker, single VCF SNP re-annotation:
     nextflow run main.nf \\
       -profile local_docker \\
       --vcf /data/sample.vcf \\
       --ref_fa /refs/hg38.fa \\
       --mode snp \\
       --hpo /data/hpo.txt \\
       --out_prefix patient_vcf_snp

  3) Local + Docker, single pre-annotated ANNOVAR SNP input:
     nextflow run main.nf \\
       -profile local_docker \\
       --annotated_snv yes \\
       --annovar_txt /data/sample.hg38_multianno.txt \\
       --vcf /data/sample.hg38_multianno.vcf \\
       --hpo /data/hpo.txt \\
       --out_prefix patient_annotated_snv

  4) Local + Docker, imported ANNOVAR SNP plus short-read SV/CNV/mito:
     nextflow run main.nf \\
       -profile local_docker \\
       --annotated_snv yes \\
       --annovar_txt /data/sample.hg38_multianno.txt \\
       --vcf /data/sample.hg38_multianno.vcf \\
       --bam /data/sample.cram \\
       --ref_fa /refs/hg38.fa \\
       --expansionhunter_variant_catalog /refs/variant_catalog.json \\
       --type short \\
       --mito yes \\
       --hpo /data/sample.hpo.txt \\
       --out_prefix patient_annotated_all_ngs

  5) CSV batch BAM with HPO inputs in note_path (PhenoTagger off):
     nextflow run main.nf \\
       -profile slurm_singularity \\
       --input_csv /data/samples.csv \\
       --bam true \\
       --ref_fa /refs/hg38.fa \\
       --expansionhunter_variant_catalog /refs/variant_catalog.json \\
       --note no

  6) Local + Singularity with custom mount sources:
     nextflow run main.nf \\
       -profile local_singularity \\
       --annovar_host_path /project/annovar \\
       --phenosv_host_path /project/train_data \\
       --bam /data/sample.bam \\
       --ref_fa /refs/hg38.fa \\
       --note /data/note.txt

NOTES
  - BAM mode requires index files (.bai for BAM, .crai for CRAM).
  - Reference FASTA index (.fai) must exist.
  - Mito mode also requires the matching `.dict` and BWA sidecars:
    `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`.
  - The DRAGEN CRAM compatibility path uses `RevertSam --RESTORE_HARDCLIPS false`.
  - For VCF single-file mode, provide --mode snp or --mode sv.
  - For single-file mode, at least one of --note <FILE> or --hpo <FILE> is required.
  - Annotated-SNV mode requires an ANNOVAR multianno TXT plus matching ANNOVAR multianno VCF.
  - Annotated-SNV plus BAM/CRAM is supported for short-read all-NGS analysis only.
  - Short-read BAM/CRAM workflows that run ExpansionHunter require
    --expansionhunter_variant_catalog <FILE>.
  
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
def raw_target_param = params.target ?: params.targeted
def clean_target = raw_target_param ? raw_target_param.trim().toLowerCase() : 'no'
def clean_note   = params.note   ? params.note.trim().toLowerCase()   : 'no'
def clean_inheritance_mode = params.inheritance_mode ? params.inheritance_mode.trim().toLowerCase() : 'ml'
def clean_include_clinvar_report = params.include_clinvar_report ? params.include_clinvar_report.trim().toLowerCase() : 'yes'
def clean_allow_unphased_comphet = params.allow_unphased_comphet ? params.allow_unphased_comphet.trim().toLowerCase() : 'no'
def clean_rankscore_softwares = params.rankscore_softwares ? params.rankscore_softwares.toString().trim() : ""
def clean_cnvnator = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : 'yes'
def clean_scramble = params.scramble ? params.scramble.toString().trim().toLowerCase() : 'no'
def clean_mito = params.mito ? params.mito.toString().trim().toLowerCase() : 'no'
def clean_annotated_snv = params.annotated_snv ? params.annotated_snv.toString().trim().toLowerCase() : 'no'

// ------------------------------------------------------------------
// 1. INPUT VALIDATION (Catching Typos)
// ------------------------------------------------------------------

// Define allowed vocabularies
def valid_modes   = ['snp', 'sv']
def valid_types   = ['ont', 'pacbio', 'short']
def valid_genomes = ['hg38', 'grch38']
def valid_inheritance_modes = ['ml', 'omim', 'gnomad']
def valid_yes_no = ['yes', 'no']
def valid_input_kinds = ['annotated_snv', 'vcf_snv', 'vcf_sv', 'bam_ngs', 'cram_ngs']
def valid_phenotype_formats = ['clinical_note', 'hpo']

def manifestHeaderColumns = []
def manifestRowsForValidation = []
def manifestUsesUnifiedSchema = false
def manifestInputKinds = [] as Set
def manifestPhenotypeFormats = [] as Set
def manifestAnnotatedRowsWithAlignment = 0
def manifestAnnotatedRowsWithoutAlignment = 0

if (params.input_csv) {
    def manifestFile = file(params.input_csv)
    if (!manifestFile.exists()) {
        error "ERROR: input CSV not found: ${params.input_csv}"
    }

    def manifestLines = manifestFile.readLines()
    if (manifestLines.isEmpty()) {
        error "ERROR: input CSV is empty: ${params.input_csv}"
    }

    manifestHeaderColumns = manifestLines[0].split(/\s*,\s*/, -1).collect { it.trim() }
    manifestUsesUnifiedSchema = manifestHeaderColumns.contains('input_kind')

    if (manifestUsesUnifiedSchema) {
        def headerIndex = [:]
        manifestHeaderColumns.eachWithIndex { column, idx -> headerIndex[column] = idx }
        def requiredHeaders = ['sample', 'input_kind', 'phenotype_path', 'phenotype_format']
        def missingHeaders = requiredHeaders.findAll { !manifestHeaderColumns.contains(it) }
        if (!missingHeaders.isEmpty()) {
            error "ERROR: Unified input CSV is missing required header(s): ${missingHeaders.join(', ')}"
        }

        def seenSamples = [] as Set
        manifestLines.drop(1).findAll { it.trim() }.eachWithIndex { line, rowIdx ->
            def values = line.split(/,/, -1).collect { it.trim() }
            if (values.size() != manifestHeaderColumns.size()) {
                error "ERROR: Unified input CSV row ${rowIdx + 2} has ${values.size()} column(s); expected ${manifestHeaderColumns.size()}."
            }

            def row = [:]
            manifestHeaderColumns.eachWithIndex { column, idx -> row[column] = values[idx] }
            def sample = row.sample
            if (!sample) {
                error "ERROR: Unified input CSV row ${rowIdx + 2} is missing a non-empty sample value."
            }
            if (seenSamples.contains(sample)) {
                error "ERROR: Unified input CSV contains duplicate sample '${sample}'. Samples must be unique in v1."
            }
            seenSamples << sample

            def inputKind = row.input_kind?.toLowerCase()
            def phenotypeFormat = row.phenotype_format?.toLowerCase()
            if (!valid_input_kinds.contains(inputKind)) {
                error "ERROR: Unified input CSV sample '${sample}' has invalid input_kind '${row.input_kind}'."
            }
            if (!valid_phenotype_formats.contains(phenotypeFormat)) {
                error "ERROR: Unified input CSV sample '${sample}' has invalid phenotype_format '${row.phenotype_format}'."
            }

            def snvTxtPath = row.containsKey('snv_txt_path') ? row.snv_txt_path : ''
            def snvVcfPath = row.containsKey('snv_vcf_path') ? row.snv_vcf_path : ''
            def vcfPath = row.containsKey('vcf_path') ? row.vcf_path : ''
            def alignmentPath = row.containsKey('alignment_path') ? row.alignment_path : ''
            def alignmentIndexPath = row.containsKey('alignment_index_path') ? row.alignment_index_path : ''

            if (inputKind == 'annotated_snv') {
                if (!snvTxtPath || !snvVcfPath) {
                    error "ERROR: Unified input CSV sample '${sample}' requires snv_txt_path and snv_vcf_path for input_kind=annotated_snv."
                }
                if (vcfPath) {
                    error "ERROR: Unified input CSV sample '${sample}' cannot mix annotated_snv fields with vcf_path."
                }
                if (alignmentPath) {
                    if (!(alignmentPath.endsWith('.bam') || alignmentPath.endsWith('.cram'))) {
                        error "ERROR: Unified input CSV sample '${sample}' requires alignment_path ending in .bam or .cram when provided with input_kind=annotated_snv."
                    }
                    manifestAnnotatedRowsWithAlignment += 1
                }
                else {
                    manifestAnnotatedRowsWithoutAlignment += 1
                }
            }
            else if (inputKind in ['vcf_snv', 'vcf_sv']) {
                if (!vcfPath) {
                    error "ERROR: Unified input CSV sample '${sample}' requires vcf_path for input_kind=${inputKind}."
                }
                if (snvTxtPath || snvVcfPath || alignmentPath || alignmentIndexPath) {
                    error "ERROR: Unified input CSV sample '${sample}' cannot mix ${inputKind} with snv_txt_path/snv_vcf_path/alignment_path/alignment_index_path."
                }
            }
            else if (inputKind in ['bam_ngs', 'cram_ngs']) {
                if (!alignmentPath) {
                    error "ERROR: Unified input CSV sample '${sample}' requires alignment_path for input_kind=${inputKind}."
                }
                if (snvTxtPath || snvVcfPath || vcfPath) {
                    error "ERROR: Unified input CSV sample '${sample}' cannot mix ${inputKind} with snv_txt_path/snv_vcf_path/vcf_path."
                }
                if (inputKind == 'bam_ngs' && !alignmentPath.endsWith('.bam')) {
                    error "ERROR: Unified input CSV sample '${sample}' requires a .bam alignment_path for input_kind=bam_ngs."
                }
                if (inputKind == 'cram_ngs' && !alignmentPath.endsWith('.cram')) {
                    error "ERROR: Unified input CSV sample '${sample}' requires a .cram alignment_path for input_kind=cram_ngs."
                }
            }

            manifestRowsForValidation << row
            manifestInputKinds << inputKind
            manifestPhenotypeFormats << phenotypeFormat
        }

        if (manifestRowsForValidation.isEmpty()) {
            error "ERROR: Unified input CSV must contain at least one sample row."
        }

        def supportedKindSets = [
            ['annotated_snv'] as Set,
            ['vcf_snv'] as Set,
            ['vcf_sv'] as Set,
            ['bam_ngs'] as Set,
            ['cram_ngs'] as Set,
            ['bam_ngs', 'cram_ngs'] as Set
        ]
        if (!supportedKindSets.any { it == manifestInputKinds }) {
            error "ERROR: Mixed unified manifests are not supported in v1. Found input kinds: ${manifestInputKinds.join(', ')}"
        }
        if (manifestAnnotatedRowsWithAlignment > 0 && manifestAnnotatedRowsWithoutAlignment > 0) {
            error "ERROR: Unified annotated_snv manifests must either provide alignment_path for every row or for none of them in v1."
        }
    }
}

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

if (!valid_yes_no.contains(clean_mito)) {
    error """
    ================================================================
    ERROR: Invalid Mito Mode
    ================================================================
    You provided: --mito "${params.mito}"

    Valid options are:
      --mito yes
      --mito no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_annotated_snv)) {
    error """
    ================================================================
    ERROR: Invalid Annotated-SNV Toggle
    ================================================================
    You provided: --annotated_snv "${params.annotated_snv}"

    Valid options are:
      --annotated_snv yes
      --annotated_snv no
    ================================================================
    """
}

// CHECK 2b: Validate CNVnator toggle
if (!valid_yes_no.contains(clean_cnvnator)) {
    error """
    ================================================================
    ERROR: Invalid CNVnator Toggle
    ================================================================
    You provided: --cnvnator "${params.cnvnator}"

    Valid options are:
      --cnvnator yes
      --cnvnator no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_scramble)) {
    error """
    ================================================================
    ERROR: Invalid SCRAMBLE Toggle
    ================================================================
    You provided: --scramble "${params.scramble}"

    Valid options are:
      --scramble yes
      --scramble no
    ================================================================
    """
}

if (clean_scramble == 'yes') {
    if (params.vcf) {
        error """
        ================================================================
        ERROR: SCRAMBLE requires BAM/CRAM input
        ================================================================
        --scramble yes is only supported for BAM/CRAM input.
        ================================================================
        """
    }

    if (clean_type != 'short') {
        error """
        ================================================================
        ERROR: SCRAMBLE is short-read only
        ================================================================
        --scramble yes currently supports only --type short.
        ================================================================
        """
    }

    if (clean_mode == 'snp') {
        error """
        ================================================================
        ERROR: SCRAMBLE is not available with --mode snp
        ================================================================
        Use --mode sv or omit --mode when running --scramble yes.
        ================================================================
        """
    }

}

if (clean_mito == 'yes') {
    if (params.vcf) {
        error """
        ================================================================
        ERROR: Mitochondrial analysis requires BAM/CRAM input
        ================================================================
        --mito yes is only supported for BAM/CRAM input in v1.
        ================================================================
        """
    }

    if (clean_type != 'short') {
        error """
        ================================================================
        ERROR: Mitochondrial analysis is short-read only
        ================================================================
        --mito yes currently supports only --type short.
        ================================================================
        """
    }

    if (clean_mode == 'sv') {
        error """
        ================================================================
        ERROR: Mitochondrial analysis is not available with --mode sv
        ================================================================
        Use --mode snp or omit --mode when running --mito yes.
        ================================================================
        """
    }
}

// CHECK 2c: Validate CNVnator bin size
if (!(params.cnvnator_bin_size.toString() ==~ /[1-9][0-9]*/)) {
    error """
    ================================================================
    ERROR: Invalid CNVnator Bin Size
    ================================================================
    You provided: --cnvnator_bin_size "${params.cnvnator_bin_size}"

    Provide a positive integer, for example:
      --cnvnator_bin_size 100
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

// CHECK 4: Validate Inheritance Assignment Mode
if (!valid_inheritance_modes.contains(clean_inheritance_mode)) {
    error """
    ================================================================
    ERROR: Invalid Inheritance Assignment Mode
    ================================================================
    You provided: --inheritance_mode "${params.inheritance_mode}"

    Valid options are:
      --inheritance_mode ml
      --inheritance_mode omim
      --inheritance_mode gnomad

    Notes:
      - ml    = ML-first with OMIM fallback (current script behavior)
      - gnomad maps to LOEUF fallback lists
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_include_clinvar_report)) {
    error """
    ================================================================
    ERROR: Invalid ClinVar Report Toggle
    ================================================================
    You provided: --include_clinvar_report "${params.include_clinvar_report}"

    Valid options are:
      --include_clinvar_report yes
      --include_clinvar_report no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_allow_unphased_comphet)) {
    error """
    ================================================================
    ERROR: Invalid Unphased CompHet Toggle
    ================================================================
    You provided: --allow_unphased_comphet "${params.allow_unphased_comphet}"

    Valid options are:
      --allow_unphased_comphet yes
      --allow_unphased_comphet no
    ================================================================
    """
}

if (clean_rankscore_softwares && clean_rankscore_softwares.split(",").every { it.trim().isEmpty() }) {
    error """
    ================================================================
    ERROR: Invalid RankScore Software List
    ================================================================
    You provided: --rankscore_softwares "${params.rankscore_softwares}"

    Provide a comma-separated software list, e.g.:
      --rankscore_softwares "REVEL,AlphaMissense,CADD_raw"
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
if (params.input_csv && !manifestUsesUnifiedSchema && !params.bam && !params.vcf) {
    error """
    ERROR: Input CSV specified, but data type is undefined.
    When using --input_csv, you must also specify the data type flag to trigger the correct logic:
      --bam true   (if CSV contains BAM/CRAM paths)
      --vcf true   (if CSV contains VCF paths)
    """
}

if (params.input_csv && manifestUsesUnifiedSchema && (params.bam || params.vcf)) {
    error """
    ERROR: Unified input CSV should not be combined with --bam true or --vcf true.
    The manifest row input_kind already selects the input path.
    """
}

// 3. Check for Missing Reference (Mandatory for BAM/CRAM)
// Logic: If BAM is being processed (either single or via CSV), Ref is required.
def unifiedManifestNeedsReference = manifestUsesUnifiedSchema && manifestRowsForValidation.any { row -> row.input_kind != 'annotated_snv' }
def annotatedHybridNeedsReference = clean_annotated_snv == 'yes' && params.bam
def annotatedManifestNeedsReference = manifestUsesUnifiedSchema && manifestRowsForValidation.any { row ->
    row.input_kind == 'annotated_snv' && row.containsKey('alignment_path') && row.alignment_path
}
def unifiedManifestNeedsExpansionHunter = manifestUsesUnifiedSchema && (
    manifestInputKinds.any { it in ['bam_ngs', 'cram_ngs'] } || annotatedManifestNeedsReference
)
def is_ref_mode = params.bam || (params.input_csv && params.bam) || (params.vcf && clean_annotated_snv != 'yes') || unifiedManifestNeedsReference || annotatedHybridNeedsReference || annotatedManifestNeedsReference
if (is_ref_mode && !params.ref_fa) {
    error """
    ERROR: Reference FASTA missing.
    This input mode requires a reference genome.
    Please specify: --ref_fa <path/to/reference.fasta>
    """
}

def is_expansionhunter_mode = clean_type == 'short' && (
    params.bam || (params.input_csv && params.bam) || unifiedManifestNeedsExpansionHunter
)
if (is_expansionhunter_mode && !params.expansionhunter_variant_catalog) {
    error """
    ERROR: ExpansionHunter variant catalog missing.
    Short-read BAM/CRAM workflows that run ExpansionHunter require:
      --expansionhunter_variant_catalog <path/to/variant_catalog.json>
    """
}

if (clean_annotated_snv == 'yes') {
    if (clean_mode == 'sv') {
        error "ERROR: Annotated-SNV mode is SNP-only. Use --mode snp or omit --mode."
    }
    if (clean_target == 'yes') {
        error "ERROR: --target yes is not supported in annotated-SNV mode."
    }
    if (!params.input_csv && !params.annovar_txt) {
        error "ERROR: Annotated-SNV single-sample mode requires --annovar_txt <FILE>."
    }
    if (!params.input_csv && !params.vcf) {
        error "ERROR: Annotated-SNV mode requires the matching --vcf <FILE>."
    }
    if (params.input_csv && !manifestUsesUnifiedSchema) {
        error "ERROR: Annotated-SNV batch mode requires the unified input CSV schema with input_kind=annotated_snv."
    }
    if (params.bam && clean_type != 'short') {
        error "ERROR: Annotated-SNV plus BAM/CRAM currently supports only --type short."
    }
}

if (manifestUsesUnifiedSchema && manifestInputKinds == (['annotated_snv'] as Set)) {
    if (manifestAnnotatedRowsWithAlignment > 0 && clean_type != 'short') {
        error "ERROR: Unified annotated_snv manifests with alignment_path currently support only --type short."
    }
    if (clean_mode == 'sv') {
        error "ERROR: Annotated-SNV manifest mode is SNP-led. Use --mode snp or omit --mode."
    }
    if (clean_target == 'yes') {
        error "ERROR: --target yes is not supported for annotated_snv manifest mode."
    }
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
    if (clean_mito == 'yes') {
        def ref_dict = file("${file(params.ref_fa).parent}/${file(params.ref_fa).baseName}.dict")
        if (!ref_dict.exists()) {
            error """
            ================================================================
            ERROR: Reference dictionary (.dict) not found
            ================================================================
            Mitochondrial analysis requires a sequence dictionary for GATK.
            Expected at: ${ref_dict}
            Create it with: gatk CreateSequenceDictionary -R ${params.ref_fa}
            ================================================================
            """
        }

        def bwa_suffixes = ['amb', 'ann', 'bwt', 'pac', 'sa']
        def missing_bwa_indexes = bwa_suffixes.findAll { suffix ->
            !file("${params.ref_fa}.${suffix}").exists()
        }
        if (!missing_bwa_indexes.isEmpty()) {
            def expected_files = missing_bwa_indexes.collect { suffix -> "${params.ref_fa}.${suffix}" }.join('\n            ')
            error """
            ================================================================
            ERROR: BWA reference sidecar files not found
            ================================================================
            Mitochondrial analysis requires a BWA-indexed reference bundle.
            Missing:
              ${expected_files}

            Create the sidecars with:
              bwa index ${params.ref_fa}
            ================================================================
            """
        }
    }
}

if (params.expansionhunter_variant_catalog) {
    def eh_catalog = file(params.expansionhunter_variant_catalog)
    if (!eh_catalog.exists()) {
        error """
        ERROR: ExpansionHunter variant catalog not found.
        Expected at: ${eh_catalog}
        Please provide a valid file with:
          --expansionhunter_variant_catalog <path/to/variant_catalog.json>
        """
    }
}

include { SINGLE_ALIGNMENT_ALL_LONGPHASE } from './subworkflows/single_alignment_all_longphase' 
include { SINGLE_ALIGNMENT_ALL_NGS } from './subworkflows/single_alignment_all_ngs'
include { SINGLE_ALIGNMENT_NGS_MITO } from './subworkflows/single_alignment_ngs_mito'
include { SINGLE_ALIGNMENT_LONG_SNP } from './subworkflows/single_alignment_long_snp'
include { SINGLE_ALIGNMENT_LONG_SV } from './subworkflows/single_alignment_long_sv'
include { SINGLE_ALIGNMENT_NGS_SNP } from './subworkflows/single_alignment_ngs_snp'
include { SINGLE_ALIGNMENT_NGS_SV } from './subworkflows/single_alignment_ngs_sv'
include { SINGLE_ALIGNMENT_VCF_SNP } from './subworkflows/single_alignment_vcf_snp'
include { SINGLE_ALIGNMENT_VCF_SV } from './subworkflows/single_alignment_vcf_sv'
include { SINGLE_ANNOTATED_SNV_SNP } from './subworkflows/single_annotated_snv_snp'
include { SINGLE_ANNOTATED_ALL_NGS } from './subworkflows/single_annotated_all_ngs'

include { INPUT_CSV_ALIGNMENT_ALL_LONGPHASE } from './subworkflows/input_csv_alignment_all_longphase'
include { INPUT_CSV_ALIGNMENT_ALL_NGS } from './subworkflows/input_csv_alignment_all_ngs'
include { INPUT_CSV_ALIGNMENT_NGS_MITO } from './subworkflows/input_csv_alignment_ngs_mito'
include { INPUT_CSV_ALIGNMENT_LONG_SNP } from './subworkflows/input_csv_alignment_long_snp'
include { INPUT_CSV_ALIGNMENT_LONG_SV } from './subworkflows/input_csv_alignment_long_sv'
include { INPUT_CSV_ALIGNMENT_NGS_SV } from './subworkflows/input_csv_alignment_ngs_sv'
include { INPUT_CSV_NGS_SNP } from './subworkflows/input_csv_alignment_ngs_snp'
include { INPUT_CSV_ALIGNMENT_VCF_SNP } from './subworkflows/input_csv_alignment_vcf_snp'
include { INPUT_CSV_ALIGNMENT_VCF_SV } from './subworkflows/input_csv_alignment_vcf_sv'
include { INPUT_CSV_ANNOTATED_SNV_SNP } from './subworkflows/input_csv_annotated_snv_snp'
include { INPUT_CSV_ANNOTATED_ALL_NGS } from './subworkflows/input_csv_annotated_all_ngs'




workflow {
	def input_age = null
	def mito_ref_fa = null
	def eh_ref_fa = null
	def eh_variant_catalog = null
	def input_annotated_snv = null
	def input_annotated_ngs = null
	def csv_manifest_mode = null
	def csv_manifest_is_note = null
	def phenotype_format = null
	if ( params.input_csv ) {
		input_age = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def out_prefix = row.sample
        def has_age_of_onset = row.containsKey('age_of_onset')
        def has_age = row.containsKey('age')
        def age_source = has_age_of_onset ? 'age_of_onset' : (has_age ? 'age' : null)
        def age_value = ''

        if (age_source != null) {
            def raw_age = row[age_source]
            age_value = raw_age == null ? '' : raw_age.toString().trim().toLowerCase()
            def age_match = (age_value =~ /^(\d+)([dmy])?$/)
            if (age_value && !age_match.matches()) {
                error """
                ERROR: Invalid age value in input CSV for sample '${out_prefix}'.
                Column '${age_source}' must be empty or match one of:
                  - <integer>
                  - <integer><unit> where unit is d/m/y
                Received: '${age_value}'
                """
            }
            if (age_value && !age_match.group(2)) {
                age_value = "${age_match.group(1)}y"
            }
        }

        return tuple(out_prefix, age_value)
    }
		if ( manifestUsesUnifiedSchema ) {
		if ( manifestInputKinds == (['annotated_snv'] as Set) ) {
		if (manifestAnnotatedRowsWithAlignment > 0) {
		input_annotated_ngs = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def bam_file = file(row.alignment_path, checkIfExists: true)
        def index_path = row.alignment_index_path ? row.alignment_index_path.toString().trim() : ''
        def bai = index_path ? file(index_path) : null
        if (bai == null || !bai.exists()) {
            def default_index = bam_file.name.endsWith('.cram')
                ? file("${bam_file}.crai")
                : file(bam_file.toString().replaceFirst(/\.bam$/, ".bai"))
            bai = default_index
            if (bam_file.name.endsWith('.cram') && !bai.exists()) {
                bai = file(bam_file.toString().replaceFirst(/\.cram$/, ".crai"))
            }
            if (bam_file.name.endsWith('.bam') && !bai.exists()) {
                bai = file("${bam_file}.bai")
            }
        }
        if (!bai.exists()) {
            error "Index file not found for ${bam_file} from unified CSV sample '${row.sample}'."
        }
        return tuple(
            row.sample,
            file(row.snv_txt_path, checkIfExists: true),
            file(row.snv_vcf_path, checkIfExists: true),
            bam_file,
            bai,
            file(row.phenotype_path, checkIfExists: true),
            row.phenotype_format.toString().trim().toLowerCase()
        )
    }
    input_annotated_snv = null
    input_vcf = null
    input_bam = null
    csv_manifest_mode = 'annotated_all_ngs'
		}
		else {
		input_annotated_snv = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        return tuple(
            row.sample,
            file(row.snv_txt_path, checkIfExists: true),
            file(row.snv_vcf_path, checkIfExists: true),
            file(row.phenotype_path, checkIfExists: true),
            row.phenotype_format.toString().trim().toLowerCase()
        )
    }
    input_annotated_ngs = null
    input_vcf = null
    input_bam = null
    csv_manifest_mode = 'annotated_snv'
		}
		}
		else if ( manifestInputKinds == (['vcf_snv'] as Set) || manifestInputKinds == (['vcf_sv'] as Set) ) {
		input_vcf = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        return tuple(
            row.sample,
            file(row.vcf_path, checkIfExists: true),
            file(row.phenotype_path, checkIfExists: true),
        )
    }
    input_bam = null
    input_annotated_snv = null
    csv_manifest_mode = manifestInputKinds.contains('vcf_sv') ? 'sv' : 'snp'
    if (manifestPhenotypeFormats.size() != 1) {
        error "ERROR: Mixed phenotype_format values are not yet supported for non-annotated CSV modes in v1."
    }
    csv_manifest_is_note = manifestPhenotypeFormats.contains('clinical_note') ? "yes" : "no"
		}
		else {
input_bam = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def bam_file = file(row.alignment_path, checkIfExists: true)
        def index_path = row.alignment_index_path ? row.alignment_index_path.toString().trim() : ''
        def bai = index_path ? file(index_path) : null
        if (bai == null || !bai.exists()) {
            def default_index = bam_file.name.endsWith('.cram')
                ? file("${bam_file}.crai")
                : file(bam_file.toString().replaceFirst(/\.bam$/, ".bai"))
            bai = default_index
            if (bam_file.name.endsWith('.cram') && !bai.exists()) {
                bai = file(bam_file.toString().replaceFirst(/\.cram$/, ".crai"))
            }
            if (bam_file.name.endsWith('.bam') && !bai.exists()) {
                bai = file("${bam_file}.bai")
            }
        }
        if (!bai.exists()) {
            error "Index file not found for ${bam_file} from unified CSV sample '${row.sample}'."
        }
        return tuple(
            row.sample,
            bam_file,
            bai,
            file(row.phenotype_path, checkIfExists: true),
        )
    }
    input_vcf = null
    input_annotated_snv = null
    csv_manifest_mode = null
    if (manifestPhenotypeFormats.size() != 1) {
        error "ERROR: Mixed phenotype_format values are not yet supported for non-annotated CSV alignment modes in v1."
    }
    csv_manifest_is_note = manifestPhenotypeFormats.contains('clinical_note') ? "yes" : "no"
		}
		}
		else if ( params.vcf ) {
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
	input_annotated_snv=null
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
            if (!bai.exists()) {
                bai = file(bam_file.toString().replaceFirst(/\.cram$/, ".crai"))
            }
        }
        
        // Ensure index actually exists before returning tuple to avoid downstream staging errors
        if (!bai.exists()) {
            error "Index file not found for ${bam_file}. Looked for BAM index paths (${bai_path}, ${bam_file}.bai) or CRAM index paths (${bam_file}.crai, ${bam_file.toString().replaceFirst(/\.cram$/, '.crai')})"
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
	input_annotated_snv=null
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
                    ? file("${file}.crai")
                    : file.parent / "${file.name}.bai"

        if (file.name.endsWith('.cram') && !index.exists()) {
            index = file(file.toString().replaceFirst(/\.cram$/, ".crai"))
        }
        if (file.name.endsWith('.bam') && !index.exists()) {
            index = file(file.toString().replaceFirst(/\.bam$/, ".bai"))
        }
        if (!index.exists()) {
            error "Index file not found for ${file}. Looked for BAM index paths (${file}.bai, ${file.toString().replaceFirst(/\.bam$/, '.bai')}) or CRAM index paths (${file}.crai, ${file.toString().replaceFirst(/\.cram$/, '.crai')})"
        }
        
        return [ file, index ]
    }
    .set { bam }
                out_prefix=Channel.value(params.out_prefix)
        }
        else if ( clean_annotated_snv == 'yes' && params.vcf != null ) {
                annovar_txt=Channel.value(file(params.annovar_txt, checkIfExists: true))
                vcf=Channel.value(file(params.vcf, checkIfExists: true))
                out_prefix=Channel.value(params.out_prefix)
        }
        else if ( params.vcf != null ) {
                vcf=Channel.value(file(params.vcf, checkIfExists: true))
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
	if (params.ref_fa != null && clean_type == 'short') {
eh_ref_fa = Channel
    .fromPath(params.ref_fa)
    .map { fa_file ->
        def fai_file = file("${fa_file}.fai")
        return [ fa_file, fai_file ]
    }
    .first()
	}
	if (params.expansionhunter_variant_catalog != null) {
eh_variant_catalog = Channel
    .fromPath(params.expansionhunter_variant_catalog)
    .first()
	}
	if (params.ref_fa != null && clean_mito == 'yes' && params.type == 'short') {
mito_ref_fa = Channel
    .fromPath(params.ref_fa)
    .map { fa_file ->
        def fai_file = file("${fa_file}.fai")
        def dict_file = file("${fa_file.parent}/${fa_file.baseName}.dict")
        def bwa_amb = file("${fa_file}.amb")
        def bwa_ann = file("${fa_file}.ann")
        def bwa_bwt = file("${fa_file}.bwt")
        def bwa_pac = file("${fa_file}.pac")
        def bwa_sa = file("${fa_file}.sa")
        return [ fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ]
    }
    .first()
	}
	// `is_note` controls whether subworkflows run phenotagger ("yes") or treat the input as HPO IDs ("no").
	def is_note = "no"
	if ( params.input_csv ) {
		if ( manifestUsesUnifiedSchema ) {
			is_note = csv_manifest_is_note
		}
		else {
			// CSV mode default: note_path is treated as clinical notes unless user explicitly sets --note no.
			is_note = (clean_note == 'no') ? "no" : "yes"
		}
	}
	else {
		// Single-file mode: `note` must be a file path; `hpo` is used directly as HPO input.
		if ( params.note != null && clean_note != 'yes' && clean_note != 'no' ) {
			note=Channel.value(params.note)
			is_note = "yes"
		}
		else if ( params.hpo != null ) {
			note=Channel.value(params.hpo)
			is_note = "no"
		}
		else {
			error """
			ERROR: Missing phenotype input for single-file mode.
			Provide one of:
			  --note <clinical_note_file>
			  --hpo  <hpo_id_file>
			"""
		}
		phenotype_format = Channel.value(is_note == "yes" ? "clinical_note" : "hpo")
	}
	output_directory_check=file(params.output_directory)
        if ( !output_directory_check.exists() ) {
	output_directory_check.mkdirs()
        }
	type=Channel.value(params.type)
	mito_contig=Channel.value(params.mito_contig)
	def inheritance_mode_script = (clean_inheritance_mode == 'gnomad') ? 'LOEUF' : 'OMIM'
	inheritance_mode=Channel.value(inheritance_mode_script)
	include_clinvar_report=Channel.value(clean_include_clinvar_report)
	allow_unphased_comphet=Channel.value(clean_allow_unphased_comphet)
	gnomad=Channel.value(params.gnomad)
        rankscore_filter=Channel.value(params.rankscore)
	rankscore_softwares=Channel.value(clean_rankscore_softwares)
	rankvar_filter=Channel.value(params.rankvar)
	phen2gene_top_n=Channel.value(params.phen2gene_filter)
	gq=Channel.value(params.gq)
	ad=Channel.value(params.ad)
	short_snp_caller = Channel.value(clean_light == 'yes' ? 'haplotypecaller' : 'deepvariant')
	long_snp_caller = Channel.value(clean_light == 'yes' ? 'nanocaller' : 'clair3')
	def target="null"
	if ( clean_target == 'yes' ) {
		target = "yes"
	}
	if ( params.input_csv ) {
        if ( input_annotated_ngs != null ) {
                INPUT_CSV_ANNOTATED_ALL_NGS(input_annotated_ngs, input_age, ref_fa, eh_ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_ref_fa, mito_contig)
        }
        else if ( input_annotated_snv != null ) {
                INPUT_CSV_ANNOTATED_SNV_SNP(input_annotated_snv, input_age, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
        }
        else if ( input_vcf != null ) {
            def batch_mode = csv_manifest_mode ?: params.mode
            if ( batch_mode == 'sv' ) {
                INPUT_CSV_ALIGNMENT_VCF_SV(input_vcf, input_age, ref_fa, phen2gene_top_n, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
            }
            else if ( batch_mode == 'snp' ) {
                INPUT_CSV_ALIGNMENT_VCF_SNP(input_vcf, input_age, ref_fa, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
            }
        }
        else if ( input_bam != null ) {
	            if ( clean_mito == 'yes' && params.type == 'short' ) {
	                INPUT_CSV_ALIGNMENT_NGS_MITO(input_bam, mito_ref_fa, mito_contig)
	            }
	            if ( params.type == 'short' ) {
	                    if ( params.mode == 'sv' ) {
	                        INPUT_CSV_ALIGNMENT_NGS_SV(input_bam, input_age, ref_fa, eh_ref_fa, eh_variant_catalog, is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else if ( params.mode == 'snp' ) {
                        INPUT_CSV_NGS_SNP(input_bam, input_age, ref_fa, eh_ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
				}
				else {
                                INPUT_CSV_ALIGNMENT_ALL_NGS(input_bam, input_age, ref_fa, eh_ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                        }
	            }
	            else { // Long reads
	                    if ( params.mode == 'sv' ) {
	                        INPUT_CSV_ALIGNMENT_LONG_SV(input_bam, input_age, ref_fa,  is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else if ( params.mode == 'snp' ) {
                        INPUT_CSV_ALIGNMENT_LONG_SNP(input_bam, input_age, ref_fa,  rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else {
                        INPUT_CSV_ALIGNMENT_ALL_LONGPHASE(input_bam, input_age, ref_fa,  rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	            }
	        }
	    }
    else { // Single File Mode
        if ( clean_annotated_snv == 'yes' && params.vcf && params.bam != null ) {
            SINGLE_ANNOTATED_ALL_NGS(bam, annovar_txt, vcf, out_prefix, ref_fa, eh_ref_fa, eh_variant_catalog, note, phenotype_format, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_ref_fa, mito_contig)
        }
        else if ( clean_annotated_snv == 'yes' && params.vcf ) {
            SINGLE_ANNOTATED_SNV_SNP(annovar_txt, vcf, out_prefix, note, phenotype_format, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
        }
        else if ( params.vcf ) {
            if ( params.mode == 'sv' ) {
                SINGLE_ALIGNMENT_VCF_SV(vcf, out_prefix, ref_fa,  note, phen2gene_top_n, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
            }
            else if ( params.mode == 'snp' ) {
                SINGLE_ALIGNMENT_VCF_SNP(vcf, out_prefix, ref_fa,  note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
            }
        }
        else if ( params.bam != null ) {
	            if ( clean_mito == 'yes' && params.type == 'short' ) {
	                SINGLE_ALIGNMENT_NGS_MITO(bam, out_prefix, mito_ref_fa, mito_contig)
	            }
	            if ( params.type == 'short' ) {
	                    if ( params.mode == 'sv' ) {
	                        SINGLE_ALIGNMENT_NGS_SV(bam, out_prefix, ref_fa, eh_ref_fa, eh_variant_catalog, note, is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else if ( params.mode == 'snp' ) {
                        SINGLE_ALIGNMENT_NGS_SNP(bam, out_prefix, ref_fa, eh_ref_fa, eh_variant_catalog, note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
			    else {
			                SINGLE_ALIGNMENT_ALL_NGS(bam, out_prefix, ref_fa, eh_ref_fa, eh_variant_catalog, note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
			    }
	            }
	            else { // Long reads
	                    if ( params.mode == 'sv' ) {
	                        SINGLE_ALIGNMENT_LONG_SV(bam, out_prefix, ref_fa,  note, is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else if ( params.mode == 'snp' ) {
                        SINGLE_ALIGNMENT_LONG_SNP(bam, out_prefix, ref_fa,  note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else {
                        SINGLE_ALIGNMENT_ALL_LONGPHASE(bam, out_prefix, ref_fa,  note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
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
