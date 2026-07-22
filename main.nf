import groovy.json.JsonSlurper
import java.security.MessageDigest

def pipelineVersion = params.pipeline_version ?: "0.3.0"

def helpMessage = """
================================================================================
  PipeVar_mito
  Version: ${pipelineVersion}
================================================================================

Rare-disease variant prioritization for short-read, long-read, and optional
mitochondrial analysis.

USAGE
  Single BAM/CRAM:
    nextflow run main.nf --bam sample.bam --ref_fa ref.fa (--note note.txt | --hpo hpo.txt) [options]

  Single VCF:
    nextflow run main.nf --vcf sample.vcf --ref_fa ref.fa --mode <snp|sv> (--note note.txt | --hpo hpo.txt) [options]

  Annotated SNV:
    nextflow run main.nf --annotated_snv yes --annovar_txt sample.hg38_multianno.txt --vcf sample.hg38_multianno.vcf (--note note.txt | --hpo hpo.txt) [options]

  Legacy CSV:
    nextflow run main.nf --input_csv samples.csv --bam true --ref_fa ref.fa [options]
    nextflow run main.nf --input_csv samples.csv --vcf true --ref_fa ref.fa --mode <snp|sv> [options]

  Unified CSV:
    nextflow run main.nf --input_csv samples.csv --ref_fa ref.fa [options]

PROFILES
  standard              SLURM + Singularity
  slurm_singularity     SLURM + Singularity
  local_singularity     Local executor + Singularity
  local_docker          Local executor + Docker

COMMON OPTIONS
  --type <ont|pacbio|short>     Sequencing type for BAM/CRAM flows
  --mode <snp|sv>               Restrict to one branch where supported
  --light <yes|no>              Use lightweight callers/models where supported
  --mito <yes|no>               Enable mitochondrial analysis for BAM/CRAM input
  --xtea <yes|no>               Enable short-read mobile-element calling
  --cnvpytor <yes|no>           Enable experimental long-read CNV calling
  --GPU <yes|no>                Enable shared GPU mode for DeepVariant GPU and PhenoGPT2
  --phenotype_extractor <STR>   phenotagger or phenogpt2 for clinical notes
  --phenogpt2_model_host_path   Complete versioned new_model directory (PhenoGPT2 notes)
  --phenogpt2_cache_host_path   Optional pre-created persistent cache directory
  --out_prefix <STRING>         Single-sample output prefix
  --output_directory <DIR>      Publish directory

NOTES
  - BAM/CRAM inputs require index files (.bai or .crai).
  - Reference FASTA index (.fai) must exist when --ref_fa is supplied.
  - CNVpytor custom references can use --cnvpytor_reference_conf <reference_genomes_conf.py>;
    built-in hg19/hg38 references are detected by CNVpytor from alignment headers.
  - Single-file mode requires one phenotype source: --note <FILE> or --hpo <FILE>.
  - PhenoGPT2 clinical-note runs require an external read-only model mount;
    HPO-only inputs do not require a PhenoGPT2 model, cache, or GPU.
  - Detailed input schemas, parameter defaults, examples, and outputs are in README.md.
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
def clean_prioritize_sv_only = params.prioritize_sv_only ? params.prioritize_sv_only.toString().trim().toLowerCase() : 'no'
def clean_common_sv_filter = params.common_sv_filter ? params.common_sv_filter.toString().trim().toLowerCase() : 'no'
def clean_denovo_filter = params.denovo_filter ? params.denovo_filter.toString().trim().toLowerCase() : 'no'
def clean_denovo_role_column = params.denovo_role_column ? params.denovo_role_column.toString().trim() : 'role'
def clean_denovo_family_column = params.denovo_family_column ? params.denovo_family_column.toString().trim() : 'family_id'
def clean_denovo_vcf_sample_column = params.denovo_vcf_sample_column ? params.denovo_vcf_sample_column.toString().trim() : 'vcf_sample'
def clean_denovo_exclude_contigs = params.denovo_exclude_contigs ? params.denovo_exclude_contigs.toString().trim() : 'MT,M,chrM,chrMT'
def clean_rankscore_softwares = params.rankscore_softwares ? params.rankscore_softwares.toString().trim() : ""
def clean_GPU = params.GPU ? params.GPU.toString().trim().toLowerCase() : 'no'
def clean_phenotype_extractor = params.phenotype_extractor ? params.phenotype_extractor.toString().trim().toLowerCase() : 'phenotagger'
def clean_phenogpt2_negation = params.phenogpt2_negation ? params.phenogpt2_negation.toString().trim().toLowerCase() : 'no'
def clean_gene_filter = params.gene ? params.gene.toString().trim() : ""
if (clean_gene_filter) {
    def gene_filter_file = file(clean_gene_filter)
    if (gene_filter_file.exists()) {
        params.gene = gene_filter_file
            .readLines()
            .collect { it.trim() }
            .findAll { it && !it.startsWith('#') }
            .join(',')
    }
}
def clean_cnvnator = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : 'yes'
def clean_cnvpytor = params.cnvpytor ? params.cnvpytor.toString().trim().toLowerCase() : 'no'
def clean_cnvpytor_baf = params.cnvpytor_baf ? params.cnvpytor_baf.toString().trim().toLowerCase() : 'yes'
def clean_cnvpytor_reference_conf = params.cnvpytor_reference_conf ? params.cnvpytor_reference_conf.toString().trim() : null
def clean_xtea = params.xtea ? params.xtea.toString().trim().toLowerCase() : 'no'
def clean_mito = params.mito ? params.mito.toString().trim().toLowerCase() : 'no'
def clean_annotated_snv = params.annotated_snv ? params.annotated_snv.toString().trim().toLowerCase() : 'no'
def clean_annotated_sv = params.annotated_sv ? params.annotated_sv.toString().trim().toLowerCase() : 'no'

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
def valid_phenotype_extractors = ['phenotagger', 'phenogpt2']

// Quote-aware parsing for manifest validation; Nextflow splitCsv performs the
// corresponding runtime parsing for channel construction.
def parseCsvRecord = { String line ->
    def values = []
    def current = new StringBuilder()
    def quoted = false
    for (int idx = 0; idx < line.size(); idx++) {
        def ch = line.charAt(idx)
        if (ch == ('"' as char)) {
            if (quoted && idx + 1 < line.size() && line.charAt(idx + 1) == ('"' as char)) {
                current.append('"')
                idx++
            }
            else {
                quoted = !quoted
            }
        }
        else if (ch == (',' as char) && !quoted) {
            values << current.toString().trim()
            current.setLength(0)
        }
        else {
            current.append(ch)
        }
    }
    if (quoted) {
        error "ERROR: unterminated quoted CSV field: ${line}"
    }
    values << current.toString().trim()
    values
}

def manifestHeaderColumns = []
def manifestRowsForValidation = []
def manifestUsesUnifiedSchema = false
def manifestInputKinds = [] as Set
def manifestPhenotypeFormats = [] as Set
def manifestAnnotatedRowsWithAlignment = 0
def manifestAnnotatedRowsWithoutAlignment = 0
def manifestAnnotatedRowsWithPreannotatedSv = 0
def manifestAnnotatedRowsWithCalledSv = 0

if (params.input_csv) {
    def manifestFile = file(params.input_csv)
    if (!manifestFile.exists()) {
        error "ERROR: input CSV not found: ${params.input_csv}"
    }

    def manifestLines = manifestFile.readLines()
    if (manifestLines.isEmpty()) {
        error "ERROR: input CSV is empty: ${params.input_csv}"
    }

    manifestHeaderColumns = parseCsvRecord(manifestLines[0])
    manifestUsesUnifiedSchema = manifestHeaderColumns.contains('input_kind')

    if (manifestUsesUnifiedSchema) {
        def requiredHeaders = ['sample', 'input_kind', 'phenotype_path', 'phenotype_format']
        def missingHeaders = requiredHeaders.findAll { !manifestHeaderColumns.contains(it) }
        if (!missingHeaders.isEmpty()) {
            error "ERROR: Unified input CSV is missing required header(s): ${missingHeaders.join(', ')}"
        }

        def seenSamples = [] as Set
        manifestLines.drop(1).findAll { it.trim() }.eachWithIndex { line, rowIdx ->
            def values = parseCsvRecord(line)
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
			if (clean_denovo_filter == 'yes' && !(sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
				error "ERROR: sample '${sample}' is not filename-safe; use only letters, numbers, '.', '_' and '-'."
			}

            def inputKind = row.input_kind?.toLowerCase()
            def phenotypeFormat = row.phenotype_format?.toLowerCase()
			def role = clean_denovo_filter == 'yes' ? (row[clean_denovo_role_column] ?: '').toString().trim().toLowerCase() : 'proband'
			def isControl = clean_denovo_filter == 'yes' && role in ['father', 'mother', 'sibling']
            if (!valid_input_kinds.contains(inputKind)) {
                error "ERROR: Unified input CSV sample '${sample}' has invalid input_kind '${row.input_kind}'."
            }
            if (!isControl && !valid_phenotype_formats.contains(phenotypeFormat)) {
                error "ERROR: Unified input CSV sample '${sample}' has invalid phenotype_format '${row.phenotype_format}'."
            }
			if (!isControl && !(row.phenotype_path ?: '').toString().trim()) {
				error "ERROR: Unified input CSV proband '${sample}' requires phenotype_path."
			}

            def snvTxtPath = row.containsKey('snv_txt_path') ? row.snv_txt_path : ''
            def snvVcfPath = row.containsKey('snv_vcf_path') ? row.snv_vcf_path : ''
            def svVcfPath = row.containsKey('sv_vcf_path') ? row.sv_vcf_path : ''
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
                    if (svVcfPath) {
                        manifestAnnotatedRowsWithPreannotatedSv += 1
                    }
                    else {
                        manifestAnnotatedRowsWithCalledSv += 1
                    }
                }
                else {
                    if (svVcfPath) {
                        error "ERROR: Unified input CSV sample '${sample}' provides sv_vcf_path without alignment_path; annotated SV import is supported only with the annotated all-NGS BAM/CRAM path."
                    }
                    manifestAnnotatedRowsWithoutAlignment += 1
                }
            }
            else if (inputKind in ['vcf_snv', 'vcf_sv']) {
                if (!vcfPath) {
                    error "ERROR: Unified input CSV sample '${sample}' requires vcf_path for input_kind=${inputKind}."
                }
                if (snvTxtPath || snvVcfPath || svVcfPath || alignmentPath || alignmentIndexPath) {
                    error "ERROR: Unified input CSV sample '${sample}' cannot mix ${inputKind} with snv_txt_path/snv_vcf_path/sv_vcf_path/alignment_path/alignment_index_path."
                }
            }
            else if (inputKind in ['bam_ngs', 'cram_ngs']) {
                if (!alignmentPath) {
                    error "ERROR: Unified input CSV sample '${sample}' requires alignment_path for input_kind=${inputKind}."
                }
                if (snvTxtPath || snvVcfPath || svVcfPath || vcfPath) {
                    error "ERROR: Unified input CSV sample '${sample}' cannot mix ${inputKind} with snv_txt_path/snv_vcf_path/sv_vcf_path/vcf_path."
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
			if (!isControl) {
				manifestPhenotypeFormats << phenotypeFormat
			}
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
        if (manifestAnnotatedRowsWithPreannotatedSv > 0 && manifestAnnotatedRowsWithCalledSv > 0) {
            error "ERROR: Unified annotated_snv manifests with alignment_path must either provide sv_vcf_path for every row or leave it blank for every row in v1."
        }
    }
}

// PhenoGPT2 resources are required only when a clinical note will actually be
// processed. HPO-only inputs continue to work without a GPU or model mount.
def singleClinicalNote = !params.input_csv && params.note != null && !(clean_note in ['yes', 'no'])
def unifiedClinicalNote = params.input_csv && manifestUsesUnifiedSchema && manifestPhenotypeFormats.contains('clinical_note')
def legacyClinicalNote = params.input_csv && !manifestUsesUnifiedSchema && !(params.note != null && clean_note == 'no')
def phenogpt2WillRun = clean_phenotype_extractor == 'phenogpt2' && (singleClinicalNote || unifiedClinicalNote || legacyClinicalNote)

def parsePositiveInteger = { value, String option ->
    try {
        def parsed = value as int
        if (parsed <= 0) {
            error "ERROR: --${option} must be a positive integer; received '${value}'."
        }
        parsed
    }
    catch (Exception ignored) {
        error "ERROR: --${option} must be a positive integer; received '${value}'."
    }
}

if (phenogpt2WillRun) {
    if (clean_GPU != 'yes') {
        error "ERROR: PhenoGPT2 will process clinical notes and requires --GPU yes."
    }
    if (clean_phenogpt2_negation != 'no') {
        error "ERROR: --phenogpt2_negation yes is not supported by the externally mounted PhenoGPT2 image."
    }

    parsePositiveInteger(params.phenogpt2_batch_size, 'phenogpt2_batch_size')
    parsePositiveInteger(params.phenogpt2_chunk_batch_size, 'phenogpt2_chunk_batch_size')
    parsePositiveInteger(params.phenogpt2_max_forks, 'phenogpt2_max_forks')
    int phenogpt2Wc
    try {
        phenogpt2Wc = params.phenogpt2_wc as int
    }
    catch (Exception ignored) {
        error "ERROR: --phenogpt2_wc must be a nonnegative integer; received '${params.phenogpt2_wc}'."
    }
    if (phenogpt2Wc < 0) {
        error "ERROR: --phenogpt2_wc must be a nonnegative integer; received '${params.phenogpt2_wc}'."
    }
    if (phenogpt2Wc > 0) {
        error "ERROR: --phenogpt2_wc > 0 is not supported until the BERT filtering model is provisioned."
    }

    def validateHostDirectory = { rawValue, String option, boolean writable ->
        if (rawValue == null || !rawValue.toString().trim()) {
            error "ERROR: --${option} is required when PhenoGPT2 processes clinical notes."
        }
        def rawPath = rawValue.toString()
        if (java.util.regex.Pattern.compile('[\\s,:\\x00-\\x1F\\x7F]').matcher(rawPath).find()) {
            error "ERROR: --${option} contains whitespace, a bind delimiter, or a control character: '${rawPath}'."
        }
        def directory = new File(rawPath)
        if (!directory.isAbsolute()) {
            error "ERROR: --${option} must be an absolute path: '${rawPath}'."
        }
        if (!directory.exists() || !directory.isDirectory()) {
            error "ERROR: --${option} must be a pre-existing directory: '${rawPath}'."
        }
        def canonical = directory.canonicalFile
        if (canonical.absolutePath != directory.absolutePath) {
            error "ERROR: --${option} must already be canonical; use '${canonical.absolutePath}'."
        }
        if (!canonical.canRead()) {
            error "ERROR: --${option} is not readable: '${canonical}'."
        }
        if (writable && !canonical.canWrite()) {
            error "ERROR: --${option} is not writable: '${canonical}'."
        }
        canonical
    }

    def modelRoot = validateHostDirectory(params.phenogpt2_model_host_path, 'phenogpt2_model_host_path', false)
    def cacheRoot = params.phenogpt2_cache_host_path ? validateHostDirectory(params.phenogpt2_cache_host_path, 'phenogpt2_cache_host_path', true) : null
    def requiredMetadata = ['config.json', 'tokenizer_config.json', 'tokenizer.json', 'chat_template.jinja', 'model.safetensors.index.json']
    def resolveModelFile = { String relativePath ->
        def relativeFile = new File(relativePath)
        if (relativeFile.isAbsolute() || relativePath.tokenize('/\\').contains('..')) {
            error "ERROR: Unsafe path in PhenoGPT2 checkpoint: '${relativePath}'."
        }
        def candidate = new File(modelRoot, relativePath)
        if (!candidate.exists()) {
            error "ERROR: Missing PhenoGPT2 checkpoint file: '${relativePath}'."
        }
        def resolved = candidate.canonicalFile
        def rootPrefix = modelRoot.absolutePath + File.separator
        if (!(resolved.absolutePath == modelRoot.absolutePath || resolved.absolutePath.startsWith(rootPrefix))) {
            error "ERROR: PhenoGPT2 checkpoint symlink escapes the model root: '${relativePath}'."
        }
        if (!resolved.isFile() || !resolved.canRead() || resolved.length() == 0L) {
            error "ERROR: PhenoGPT2 checkpoint file is empty, unreadable, or not regular: '${relativePath}'."
        }
        resolved
    }

    def metadataFiles = requiredMetadata.collect { resolveModelFile(it) }
    def checkpointIndex
    try {
        checkpointIndex = new JsonSlurper().parse(metadataFiles[-1])
    }
    catch (Exception exc) {
        error "ERROR: Invalid PhenoGPT2 checkpoint index: ${exc.message}"
    }
    if (!(checkpointIndex.weight_map instanceof Map) || checkpointIndex.weight_map.isEmpty()) {
        error "ERROR: PhenoGPT2 checkpoint index has no nonempty weight_map."
    }
    if (checkpointIndex.weight_map.values().any { !(it instanceof String) || !it }) {
        error "ERROR: PhenoGPT2 checkpoint index contains an invalid shard path."
    }
    def shardNames = checkpointIndex.weight_map.values().unique().sort()
    def shardFiles = shardNames.collect { resolveModelFile(it) }

    def digest = MessageDigest.getInstance('SHA-256')
    metadataFiles.each { metadata ->
        digest.update(metadata.name.getBytes('UTF-8'))
        digest.update(metadata.bytes)
    }
    shardNames.eachWithIndex { shardName, idx ->
        def shard = shardFiles[idx]
        digest.update(shardName.getBytes('UTF-8'))
        digest.update(shard.length().toString().getBytes('UTF-8'))
        digest.update(shard.lastModified().toString().getBytes('UTF-8'))
    }
    params.phenogpt2_model_host_path = modelRoot.absolutePath
    params.phenogpt2_cache_host_path = cacheRoot?.absolutePath
    params.phenogpt2_model_fingerprint = digest.digest().collect { String.format('%02x', it & 0xff) }.join()
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

if (!valid_yes_no.contains(clean_GPU)) {
    error """
    ================================================================
    ERROR: Invalid GPU Toggle
    ================================================================
    You provided: --GPU "${params.GPU}"

    Valid options are:
      --GPU yes
      --GPU no
    ================================================================
    """
}

if (!valid_phenotype_extractors.contains(clean_phenotype_extractor)) {
    error """
    ================================================================
    ERROR: Invalid Phenotype Extractor
    ================================================================
    You provided: --phenotype_extractor "${params.phenotype_extractor}"

    Valid options are:
      --phenotype_extractor phenotagger
      --phenotype_extractor phenogpt2
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

if (!valid_yes_no.contains(clean_annotated_sv)) {
    error """
    ================================================================
    ERROR: Invalid Annotated-SV Toggle
    ================================================================
    You provided: --annotated_sv "${params.annotated_sv}"

    Valid options are:
      --annotated_sv yes
      --annotated_sv no
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

if (!valid_yes_no.contains(clean_xtea)) {
    error """
    ================================================================
    ERROR: Invalid xTEA Toggle
    ================================================================
    You provided: --xtea "${params.xtea}"

    Valid options are:
      --xtea yes
      --xtea no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_cnvpytor)) {
    error """
    ================================================================
    ERROR: Invalid CNVpytor Toggle
    ================================================================
    You provided: --cnvpytor "${params.cnvpytor}"

    Valid options are:
      --cnvpytor yes
      --cnvpytor no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_cnvpytor_baf)) {
    error """
    ================================================================
    ERROR: Invalid CNVpytor BAF Toggle
    ================================================================
    You provided: --cnvpytor_baf "${params.cnvpytor_baf}"

    Valid options are:
      --cnvpytor_baf yes
      --cnvpytor_baf no
    ================================================================
    """
}

if (clean_xtea == 'yes') {
    if (params.vcf) {
        error """
        ================================================================
        ERROR: xTEA requires BAM/CRAM input
        ================================================================
        --xtea yes is only supported for BAM/CRAM input.
        ================================================================
        """
    }

    if (clean_type != 'short') {
        error """
        ================================================================
        ERROR: xTEA is short-read only
        ================================================================
        --xtea yes currently supports only --type short.
        ================================================================
        """
    }

    if (clean_mode == 'snp') {
        error """
        ================================================================
        ERROR: xTEA is not available with --mode snp
        ================================================================
        Use --mode sv or omit --mode when running --xtea yes.
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

    if (clean_mode == 'sv') {
        error """
        ================================================================
        ERROR: Mitochondrial analysis is not available with --mode sv
        ================================================================
        Use --mode snp or omit --mode when running --mito yes.
        ================================================================
        """
    }

    if (clean_type != 'short' && clean_light == 'yes') {
        error """
        ================================================================
        ERROR: Long-read mitochondrial analysis requires Clair3
        ================================================================
        --mito yes is not available when --light yes selects NanoCaller.
        Use standard long-read mode or --type short.
        ================================================================
        """
    }
}

if (clean_cnvpytor == 'yes') {
    if (params.vcf) {
        error """
        ================================================================
        ERROR: CNVpytor requires BAM/CRAM input
        ================================================================
        --cnvpytor yes is only supported for BAM/CRAM input.
        ================================================================
        """
    }

    if (clean_type == 'short') {
        error """
        ================================================================
        ERROR: CNVpytor is long-read only
        ================================================================
        --cnvpytor yes currently supports only --type ont or --type pacbio.
        ================================================================
        """
    }

    if (clean_mode == 'snp') {
        error """
        ================================================================
        ERROR: CNVpytor is not available with --mode snp
        ================================================================
        Use --mode sv or omit --mode when running --cnvpytor yes.
        ================================================================
        """
    }

    if (clean_mode == 'sv' && clean_cnvpytor_baf == 'yes') {
        println "WARNING: --cnvpytor_baf yes was requested with --mode sv; CNVpytor will run in RD-only mode for the SV-only long-read branch."
    }

    if (clean_cnvpytor_reference_conf && !file(clean_cnvpytor_reference_conf).exists()) {
        error """
        ================================================================
        ERROR: CNVpytor Reference Config Not Found
        ================================================================
        You provided: --cnvpytor_reference_conf "${params.cnvpytor_reference_conf}"

        Provide an existing CNVpytor reference_genomes_conf.py file, or omit
        --cnvpytor_reference_conf when using built-in references such as hg19
        or hg38.
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

if (!(params.cnvpytor_primary_bin.toString() ==~ /[1-9][0-9]*/)) {
    error """
    ================================================================
    ERROR: Invalid CNVpytor Primary Bin
    ================================================================
    You provided: --cnvpytor_primary_bin "${params.cnvpytor_primary_bin}"

    Provide a positive integer, for example:
      --cnvpytor_primary_bin 100000
    ================================================================
    """
}

if (!(params.cnvpytor_min_size.toString() ==~ /[1-9][0-9]*/)) {
    error """
    ================================================================
    ERROR: Invalid CNVpytor Minimum Size
    ================================================================
    You provided: --cnvpytor_min_size "${params.cnvpytor_min_size}"

    Provide a positive integer, for example:
      --cnvpytor_min_size 100000
    ================================================================
    """
}

def clean_cnvpytor_bin_sizes = params.cnvpytor_bin_sizes ? params.cnvpytor_bin_sizes.toString().trim() : "100000"
if (!clean_cnvpytor_bin_sizes || !(clean_cnvpytor_bin_sizes ==~ /[1-9][0-9]*(\s+[1-9][0-9]*)*/)) {
    error """
    ================================================================
    ERROR: Invalid CNVpytor Bin Sizes
    ================================================================
    You provided: --cnvpytor_bin_sizes "${params.cnvpytor_bin_sizes}"

    Provide one or more positive integers separated by spaces, for example:
      --cnvpytor_bin_sizes "100000"
      --cnvpytor_bin_sizes "50000 100000"
    ================================================================
    """
}

params.mode = clean_mode
params.type = clean_type
params.light = clean_light
params.genome = clean_genome
params.cnvnator = clean_cnvnator
params.cnvpytor = clean_cnvpytor
params.cnvpytor_baf = clean_cnvpytor_baf
params.cnvpytor_bin_sizes = clean_cnvpytor_bin_sizes
params.cnvpytor_reference_conf = clean_cnvpytor_reference_conf
params.xtea = clean_xtea
params.mito = clean_mito
params.annotated_snv = clean_annotated_snv
params.annotated_sv = clean_annotated_sv

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

if (!valid_yes_no.contains(clean_prioritize_sv_only)) {
    error """
    ================================================================
    ERROR: Invalid SV-only final prioritization toggle
    ================================================================
    You provided: --prioritize_sv_only "${params.prioritize_sv_only}"

    Valid options are:
      --prioritize_sv_only yes
      --prioritize_sv_only no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_common_sv_filter)) {
    error """
    ================================================================
    ERROR: Invalid Common SV Filter Toggle
    ================================================================
    You provided: --common_sv_filter "${params.common_sv_filter}"

    Valid options are:
      --common_sv_filter yes
      --common_sv_filter no
    ================================================================
    """
}

if (!valid_yes_no.contains(clean_denovo_filter)) {
    error "ERROR: Invalid --denovo_filter '${params.denovo_filter}'. Use 'yes' or 'no'."
}

[clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column].each { columnName ->
    if (!(columnName ==~ /[A-Za-z_][A-Za-z0-9_.-]*/)) {
        error "ERROR: invalid de novo CSV column name '${columnName}'."
    }
}

if (!(params.denovo_sv_min_reciprocal_overlap.toString() ==~ /0\.[0-9]*[1-9][0-9]*|1(\.0+)?/)) {
    error "ERROR: Invalid --denovo_sv_min_reciprocal_overlap '${params.denovo_sv_min_reciprocal_overlap}'. Provide a value greater than 0 and at most 1."
}

if (!(params.phenosv_score.toString() ==~ /([0-9]+([.][0-9]+)?|[.][0-9]+)/)) {
    error """
    ================================================================
    ERROR: Invalid PhenoSV Score Threshold
    ================================================================
    You provided: --phenosv_score "${params.phenosv_score}"

    Provide a numeric threshold, for example:
      --phenosv_score 0.50
    ================================================================
    """
}

if (!(params.common_sv_af.toString() ==~ /([0-9]+([.][0-9]+)?|[.][0-9]+)/)) {
    error "ERROR: Invalid --common_sv_af '${params.common_sv_af}'. Provide a numeric threshold, for example 0.01."
}

if (!(params.common_sv_reciprocal_overlap.toString() ==~ /([0-9]+([.][0-9]+)?|[.][0-9]+)/)) {
    error "ERROR: Invalid --common_sv_reciprocal_overlap '${params.common_sv_reciprocal_overlap}'. Provide a numeric threshold, for example 0.5."
}

if (!(params.common_sv_distance.toString() ==~ /[0-9]+/)) {
    error "ERROR: Invalid --common_sv_distance '${params.common_sv_distance}'. Provide an integer distance, for example 1000."
}

if (!(params.common_sv_ins_distance.toString() ==~ /[0-9]+/)) {
    error "ERROR: Invalid --common_sv_ins_distance '${params.common_sv_ins_distance}'. Provide an integer distance, for example 500."
}

if (!(params.common_sv_ins_identity.toString() ==~ /([0-9]+([.][0-9]+)?|[.][0-9]+)/)) {
    error "ERROR: Invalid --common_sv_ins_identity '${params.common_sv_ins_identity}'. Provide a numeric threshold, for example 0.8."
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

if (clean_denovo_filter == 'yes' && !params.input_csv) {
    error "ERROR: --denovo_filter yes is supported only with --input_csv."
}

if (clean_denovo_filter == 'yes') {
    def pedigreeFile = file(params.input_csv)
    def pedigreeLines = pedigreeFile.readLines().findAll { it.trim() }
    if (pedigreeLines.isEmpty()) {
        error "ERROR: de novo pedigree CSV is empty: ${params.input_csv}"
    }
    def pedigreeHeader = parseCsvRecord(pedigreeLines[0])
    def sampleIndex = pedigreeHeader.indexOf('sample')
    def roleIndex = pedigreeHeader.indexOf(clean_denovo_role_column)
    def familyIndex = pedigreeHeader.indexOf(clean_denovo_family_column)
    if (sampleIndex < 0 || roleIndex < 0) {
        error "ERROR: --denovo_filter yes requires CSV columns 'sample' and '${clean_denovo_role_column}'."
    }

    def allowedRoles = ['proband', 'father', 'mother', 'sibling'] as Set
    def pedigreeFamilies = [:]
	def pedigreeSamples = [] as Set
    pedigreeLines.drop(1).eachWithIndex { line, rowIndex ->
        def values = parseCsvRecord(line)
        if (values.size() != pedigreeHeader.size()) {
            error "ERROR: de novo CSV row ${rowIndex + 2} has ${values.size()} column(s); expected ${pedigreeHeader.size()}."
        }
        def sample = values[sampleIndex]
        def role = values[roleIndex].toLowerCase()
        if (!sample || !allowedRoles.contains(role)) {
            error "ERROR: invalid de novo pedigree row ${rowIndex + 2}; sample must be non-empty and role must be father, mother, proband, or sibling."
        }
		if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
			error "ERROR: sample '${sample}' is not filename-safe; use only letters, numbers, '.', '_' and '-'."
		}
		if (pedigreeSamples.contains(sample)) {
			error "ERROR: duplicate sample '${sample}' in de novo CSV. Sample keys must be globally unique."
		}
		pedigreeSamples << sample
        if (role != 'sibling') {
            def familyId = familyIndex >= 0 && values[familyIndex] ? values[familyIndex] : 'family_1'
            if (!pedigreeFamilies.containsKey(familyId)) {
                pedigreeFamilies[familyId] = [:]
            }
            def family = pedigreeFamilies[familyId]
            if (family.containsKey(role)) {
                error "ERROR: family '${familyId}' has duplicate role '${role}' (${family[role]}, ${sample})."
            }
            family[role] = sample
        }
    }
    if (pedigreeFamilies.isEmpty()) {
        error "ERROR: --denovo_filter yes requires at least one proband and one parent."
    }
    pedigreeFamilies.each { familyId, family ->
        if (!family.containsKey('proband')) {
            error "ERROR: family '${familyId}' is missing a proband."
        }
        if (!family.containsKey('father') && !family.containsKey('mother')) {
            error "ERROR: family '${familyId}' must include father and/or mother."
        }
    }
	def requestedSexColumn = params.sex_column ? params.sex_column.toString().trim() : 'sex'
	if (requestedSexColumn != 'sex' && !pedigreeHeader.contains(requestedSexColumn)) {
		error "ERROR: requested --sex_column '${requestedSexColumn}' is absent from the input CSV."
	}
}

if (clean_annotated_snv == 'yes') {
    if (clean_mode == 'sv') {
        error "ERROR: Annotated-SNV mode is SNP-led. Use --mode snp or omit --mode."
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

if (clean_annotated_sv == 'yes') {
    if (clean_annotated_snv != 'yes') {
        error "ERROR: --annotated_sv yes is currently supported only with --annotated_snv yes."
    }
    if (!params.input_csv && !params.annovar_sv_vcf) {
        error "ERROR: Annotated-SV single-sample mode requires --annovar_sv_vcf <FILE>."
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
def is_bam_mode = params.bam || (params.input_csv && params.bam) || unifiedManifestNeedsReference || annotatedHybridNeedsReference || annotatedManifestNeedsReference
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
    if (clean_mito == 'yes' && clean_type == 'short') {
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

def is_expansionhunter_mode = clean_type == 'short' && (
    params.bam || (params.input_csv && params.bam) || unifiedManifestNeedsExpansionHunter
)
def default_eh_catalog = clean_genome == 'grch38' ? "${projectDir}/data/variant_catalog_grch38.json" : "${projectDir}/data/variant_catalog.json"
def selected_eh_catalog = params.expansionhunter_variant_catalog ?: default_eh_catalog

if (is_expansionhunter_mode) {
    def eh_catalog = file(selected_eh_catalog)
    if (!eh_catalog.exists()) {
        error """
        ERROR: ExpansionHunter variant catalog not found.
        Expected at: ${eh_catalog}
        Provide a valid file with:
          --expansionhunter_variant_catalog <path/to/variant_catalog.json>
        or ensure the default PipeVar_mito/data catalog exists for --genome ${clean_genome}.
        """
    }
}

include { SINGLE_ALIGNMENT_ALL_LONGPHASE } from './subworkflows/single_alignment_all_longphase' 
include { SINGLE_ALIGNMENT_ALL_NGS } from './subworkflows/single_alignment_all_ngs'
include { SINGLE_ALIGNMENT_LONG_MITO } from './subworkflows/single_alignment_long_mito'
include { SINGLE_ALIGNMENT_NGS_MITO } from './subworkflows/single_alignment_ngs_mito'
include { SINGLE_ALIGNMENT_LONG_SNP } from './subworkflows/single_alignment_long_snp'
include { SINGLE_ALIGNMENT_LONG_SV } from './subworkflows/single_alignment_long_sv'
include { SINGLE_ALIGNMENT_NGS_SNP } from './subworkflows/single_alignment_ngs_snp'
include { SINGLE_ALIGNMENT_NGS_SV } from './subworkflows/single_alignment_ngs_sv'
include { SINGLE_ALIGNMENT_VCF_SNP } from './subworkflows/single_alignment_vcf_snp'
include { SINGLE_ALIGNMENT_VCF_SV } from './subworkflows/single_alignment_vcf_sv'
include { SINGLE_ANNOTATED_SNV_SNP } from './subworkflows/single_annotated_snv_snp'
include { SINGLE_ANNOTATED_ALL_NGS } from './subworkflows/single_annotated_all_ngs'
include { SINGLE_ANNOTATED_SNV_CALLED_SV_NGS } from './subworkflows/single_annotated_snv_called_sv_ngs'

include { INPUT_CSV_ALIGNMENT_ALL_LONGPHASE } from './subworkflows/input_csv_alignment_all_longphase'
include { INPUT_CSV_ALIGNMENT_ALL_NGS } from './subworkflows/input_csv_alignment_all_ngs'
include { INPUT_CSV_ALIGNMENT_LONG_MITO } from './subworkflows/input_csv_alignment_long_mito'
include { INPUT_CSV_ALIGNMENT_NGS_MITO } from './subworkflows/input_csv_alignment_ngs_mito'
include { INPUT_CSV_ALIGNMENT_LONG_SNP } from './subworkflows/input_csv_alignment_long_snp'
include { INPUT_CSV_ALIGNMENT_LONG_SV } from './subworkflows/input_csv_alignment_long_sv'
include { INPUT_CSV_ALIGNMENT_NGS_SV } from './subworkflows/input_csv_alignment_ngs_sv'
include { INPUT_CSV_NGS_SNP } from './subworkflows/input_csv_alignment_ngs_snp'
include { INPUT_CSV_ALIGNMENT_VCF_SNP } from './subworkflows/input_csv_alignment_vcf_snp'
include { INPUT_CSV_ALIGNMENT_VCF_SV } from './subworkflows/input_csv_alignment_vcf_sv'
include { INPUT_CSV_ANNOTATED_SNV_SNP } from './subworkflows/input_csv_annotated_snv_snp'
include { INPUT_CSV_ANNOTATED_ALL_NGS } from './subworkflows/input_csv_annotated_all_ngs'
include { INPUT_CSV_ANNOTATED_SNV_CALLED_SV_NGS } from './subworkflows/input_csv_annotated_snv_called_sv_ngs'




workflow {
	def input_meta = null
	def mito_ref_fa = null
	def eh_variant_catalog = null
	def annovar_txt = null
	def annovar_sv_vcf = null
	def input_annotated_snv = null
	def input_annotated_ngs = null
	def input_annotated_called_ngs = null
	def csv_manifest_mode = null
	def csv_manifest_is_note = null
	def denovo_pedigree = Channel.value("null")
	def phenotypeFileForRow = { row, columnName ->
		def role = clean_denovo_filter == 'yes' ? (row[clean_denovo_role_column] ?: '').toString().trim().toLowerCase() : 'proband'
		if (clean_denovo_filter == 'yes' && role != 'proband') {
			return 'null'
		}
		def rawPath = (row[columnName] ?: '').toString().trim()
		if (!rawPath) {
			error "ERROR: proband '${row.sample}' requires a non-empty ${columnName}."
		}
		file(rawPath, checkIfExists: true)
	}
	if ( params.input_csv ) {
		if ( clean_denovo_filter == 'yes' ) {
			denovo_pedigree = Channel.fromPath(params.input_csv).first()
		}
		input_meta = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
	.filter { row ->
		if (clean_denovo_filter != 'yes') {
			return true
		}
		return (row[clean_denovo_role_column] ?: '').toString().trim().toLowerCase() == 'proband'
	}
    .map { row ->
        def out_prefix = row.sample
        def has_age_of_onset = row.containsKey('age_of_onset')
        def has_age = row.containsKey('age')
        def age_source = has_age_of_onset ? 'age_of_onset' : (has_age ? 'age' : null)
        def age_value = ''
        def sex_column = params.sex_column ? params.sex_column.toString() : 'sex'
        def sex_value = 'unknown'

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

        if (row.containsKey(sex_column)) {
            def raw_sex = row[sex_column]
            sex_value = raw_sex == null ? 'unknown' : raw_sex.toString().trim().toLowerCase()
            if (!sex_value) {
                sex_value = 'unknown'
            }
            if (!(sex_value in ['unknown', 'male', 'female'])) {
                error """
                ERROR: Invalid sex value in input CSV for sample '${out_prefix}'.
                Column '${sex_column}' must be empty, 'unknown', 'male', or 'female'.
                Received: '${raw_sex}'
                """
            }
        }

        return tuple(out_prefix, age_value, sex_value)
    }
		if ( manifestUsesUnifiedSchema ) {
		if ( manifestInputKinds == (['annotated_snv'] as Set) ) {
		if (manifestAnnotatedRowsWithPreannotatedSv > 0) {
		input_annotated_ngs = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def bam_file = file(row.alignment_path, checkIfExists: true)
        def index_path = row.alignment_index_path ? row.alignment_index_path.toString().trim() : ''
        def bai = index_path ? file(index_path) : null
        if (bai == null || !bai.exists()) {
            bai = bam_file.name.endsWith('.cram') ? file("${bam_file}.crai") : file(bam_file.toString().replaceFirst(/\.bam$/, ".bai"))
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
            file(row.sv_vcf_path, checkIfExists: true),
            bam_file,
            bai,
            phenotypeFileForRow(row, 'phenotype_path'),
            (row.phenotype_format ?: '').toString().trim().toLowerCase()
        )
    }
    input_annotated_snv = null
    input_annotated_called_ngs = null
    input_vcf = null
    input_bam = null
    csv_manifest_mode = 'annotated_all_ngs'
		}
		else if (manifestAnnotatedRowsWithCalledSv > 0) {
		input_annotated_called_ngs = Channel
    .fromPath( params.input_csv )
    .splitCsv( header:true )
    .map { row ->
        def bam_file = file(row.alignment_path, checkIfExists: true)
        def index_path = row.alignment_index_path ? row.alignment_index_path.toString().trim() : ''
        def bai = index_path ? file(index_path) : null
        if (bai == null || !bai.exists()) {
            bai = bam_file.name.endsWith('.cram') ? file("${bam_file}.crai") : file(bam_file.toString().replaceFirst(/\.bam$/, ".bai"))
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
            phenotypeFileForRow(row, 'phenotype_path'),
            (row.phenotype_format ?: '').toString().trim().toLowerCase()
        )
    }
    input_annotated_snv = null
    input_annotated_ngs = null
    input_vcf = null
    input_bam = null
    csv_manifest_mode = 'annotated_snv_called_sv_ngs'
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
            phenotypeFileForRow(row, 'phenotype_path'),
            (row.phenotype_format ?: '').toString().trim().toLowerCase()
        )
    }
    input_annotated_ngs = null
    input_annotated_called_ngs = null
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
            phenotypeFileForRow(row, 'phenotype_path'),
        )
    }
    input_bam = null
    input_annotated_snv = null
    input_annotated_ngs = null
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
            bai = bam_file.name.endsWith('.cram') ? file("${bam_file}.crai") : file(bam_file.toString().replaceFirst(/\.bam$/, ".bai"))
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
            phenotypeFileForRow(row, 'phenotype_path'),
        )
    }
    input_vcf = null
    input_annotated_snv = null
    input_annotated_ngs = null
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
        def note_file = phenotypeFileForRow(row, 'note_path')
	def out_prefix = row.sample

        // Return a List/tuple in a specific order
        return tuple(
            row.sample,
            vcf_file,
            note_file,
        )
    			}
	input_bam=null
	input_annotated_snv=null
	input_annotated_ngs=null
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

        def note_file = phenotypeFileForRow(row, 'note_path')
	def out_prefix = row.sample

        // Return the tuple
        return tuple(
            out_prefix,
            bam_file,
            bai,
            note_file,
        )
                        }
	input_vcf=null
	input_annotated_snv=null
	input_annotated_ngs=null

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
        else if ( params.vcf != null ) {
                vcf=Channel.value(params.vcf)
                out_prefix=Channel.value(params.out_prefix)
        }
	if (clean_annotated_snv == 'yes' && !params.input_csv) {
		annovar_txt = Channel.value(params.annovar_txt)
		vcf = Channel.value(params.vcf)
		if (clean_annotated_sv == 'yes') {
			annovar_sv_vcf = Channel.value(params.annovar_sv_vcf)
		}
	}
	//Need to add dict file in case there is we are running haplotypecaller!!!!
	if (params.ref_fa != null && clean_light == 'yes' && clean_type == 'short' ) {
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
	if (params.ref_fa != null && clean_mito == 'yes' && clean_type == 'short') {
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
	if (is_expansionhunter_mode) {
eh_variant_catalog = Channel
    .fromPath(selected_eh_catalog)
    .first()
	}
	// `is_note` controls whether subworkflows run phenotagger ("yes") or treat the input as HPO IDs ("no").
	def is_note = "no"
	if ( params.input_csv ) {
		// CSV mode default: note_path is treated as clinical notes unless user explicitly sets --note no.
		is_note = csv_manifest_is_note != null ? csv_manifest_is_note : ((params.note != null && clean_note == 'no') ? "no" : "yes")
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
	}
	output_directory_check=file(params.output_directory)
        if ( !output_directory_check.exists() ) {
	output_directory_check.mkdirs()
        }
	type=Channel.value(clean_type)
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
	mito_input_bam = input_bam
	if (clean_denovo_filter == 'yes' && input_bam != null) {
		proband_meta_keys = input_meta.map { sample, age_of_onset, sex -> tuple(sample, true) }
		mito_input_bam = input_bam.join(proband_meta_keys, failOnDuplicate: true).map { sample, bam_file, bai_file, phenotype_file, marker ->
			tuple(sample, bam_file, bai_file, phenotype_file)
		}
	}
	if ( params.input_csv ) {
        if ( input_annotated_ngs != null ) {
            INPUT_CSV_ANNOTATED_ALL_NGS(input_annotated_ngs, input_meta, ref_fa, ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_ref_fa, mito_contig, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
        }
        else if ( input_annotated_called_ngs != null ) {
            INPUT_CSV_ANNOTATED_SNV_CALLED_SV_NGS(input_annotated_called_ngs, input_meta, ref_fa, ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_ref_fa, mito_contig, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
        }
        else if ( input_annotated_snv != null ) {
            INPUT_CSV_ANNOTATED_SNV_SNP(input_annotated_snv, input_meta, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
        }
        else if ( input_vcf != null ) {
            def effective_csv_mode = clean_mode ?: csv_manifest_mode
            if ( effective_csv_mode == 'sv' ) {
                INPUT_CSV_ALIGNMENT_VCF_SV(input_vcf, input_meta, ref_fa, phen2gene_top_n, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
            }
            else if ( effective_csv_mode == 'snp' ) {
                INPUT_CSV_ALIGNMENT_VCF_SNP(input_vcf, input_meta, ref_fa, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
            }
        }
	        else if ( input_bam != null ) {
		            if ( clean_type == 'short' ) {
		                    mito_report_tsv = Channel.empty()
		                    if ( clean_mito == 'yes' ) {
		                        INPUT_CSV_ALIGNMENT_NGS_MITO(mito_input_bam, mito_ref_fa, mito_contig)
		                        mito_report_tsv = INPUT_CSV_ALIGNMENT_NGS_MITO.out.prioritized_tsv
		                    }
		                    if ( clean_mode == 'sv' ) {
		                        INPUT_CSV_ALIGNMENT_NGS_SV(input_bam, input_meta, ref_fa, eh_variant_catalog, is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
		                    }
	                    else if ( clean_mode == 'snp' ) {
                        INPUT_CSV_NGS_SNP(input_bam, input_meta, ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
				}
				else {
                                INPUT_CSV_ALIGNMENT_ALL_NGS(input_bam, input_meta, ref_fa, eh_variant_catalog, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_report_tsv, clean_mito, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
	                        }
	            }
	            else { // Long reads
	                    mito_report_tsv = Channel.empty()
	                    if ( clean_mito == 'yes' ) {
	                        INPUT_CSV_ALIGNMENT_LONG_MITO(mito_input_bam, ref_fa, mito_contig)
	                        mito_report_tsv = INPUT_CSV_ALIGNMENT_LONG_MITO.out.prioritized_tsv
	                    }
	                    if ( clean_mode == 'sv' ) {
	                        INPUT_CSV_ALIGNMENT_LONG_SV(input_bam, input_meta, ref_fa,  is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
	                    }
	                    else if ( clean_mode == 'snp' ) {
                        INPUT_CSV_ALIGNMENT_LONG_SNP(input_bam, input_meta, ref_fa,  rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
	                    }
	                    else {
                        INPUT_CSV_ALIGNMENT_ALL_LONGPHASE(input_bam, input_meta, ref_fa,  rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_report_tsv, clean_mito, clean_denovo_filter, denovo_pedigree, clean_denovo_role_column, clean_denovo_family_column, clean_denovo_vcf_sample_column, clean_denovo_exclude_contigs, params.denovo_sv_min_reciprocal_overlap)
	                    }
	            }
	        }
	    }
    else { // Single File Mode
        if ( clean_annotated_snv == 'yes' && params.vcf && params.bam != null ) {
            if ( clean_annotated_sv == 'yes' ) {
                SINGLE_ANNOTATED_ALL_NGS(bam, annovar_txt, vcf, annovar_sv_vcf, out_prefix, ref_fa, ref_fa, eh_variant_catalog, note, is_note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_ref_fa, mito_contig)
            }
            else {
                SINGLE_ANNOTATED_SNV_CALLED_SV_NGS(bam, annovar_txt, vcf, out_prefix, ref_fa, ref_fa, eh_variant_catalog, note, is_note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_ref_fa, mito_contig)
            }
        }
        else if ( clean_annotated_snv == 'yes' && params.vcf ) {
            SINGLE_ANNOTATED_SNV_SNP(annovar_txt, vcf, out_prefix, note, is_note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
        }
        else if ( params.vcf ) {
            if ( clean_mode == 'sv' ) {
                SINGLE_ALIGNMENT_VCF_SV(vcf, out_prefix, ref_fa,  note, phen2gene_top_n, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
            }
            else if ( clean_mode == 'snp' ) {
                SINGLE_ALIGNMENT_VCF_SNP(vcf, out_prefix, ref_fa,  note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
            }
        }
        else if ( params.bam != null ) {
	            if ( clean_type == 'short' ) {
	                    mito_report_tsv = Channel.empty()
	                    if ( clean_mito == 'yes' ) {
	                        SINGLE_ALIGNMENT_NGS_MITO(bam, out_prefix, mito_ref_fa, mito_contig)
	                        mito_report_tsv = SINGLE_ALIGNMENT_NGS_MITO.out.prioritized_tsv
	                    }
	                    if ( clean_mode == 'sv' ) {
	                        SINGLE_ALIGNMENT_NGS_SV(bam, out_prefix, ref_fa, eh_variant_catalog, note, is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else if ( clean_mode == 'snp' ) {
                        SINGLE_ALIGNMENT_NGS_SNP(bam, out_prefix, ref_fa, eh_variant_catalog, note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
			    else {
			                SINGLE_ALIGNMENT_ALL_NGS(bam, out_prefix, ref_fa, eh_variant_catalog, note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, short_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_report_tsv, clean_mito)
			    }
	            }
	            else { // Long reads
	                    mito_report_tsv = Channel.empty()
	                    if ( clean_mito == 'yes' ) {
	                        SINGLE_ALIGNMENT_LONG_MITO(bam, out_prefix, ref_fa, mito_contig)
	                        mito_report_tsv = SINGLE_ALIGNMENT_LONG_MITO.out.prioritized_tsv
	                    }
	                    if ( clean_mode == 'sv' ) {
	                        SINGLE_ALIGNMENT_LONG_SV(bam, out_prefix, ref_fa,  note, is_note, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else if ( clean_mode == 'snp' ) {
                        SINGLE_ALIGNMENT_LONG_SNP(bam, out_prefix, ref_fa,  note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	                    }
	                    else {
                        SINGLE_ALIGNMENT_ALL_LONGPHASE(bam, out_prefix, ref_fa,  note, rankscore_filter, rankscore_softwares, phen2gene_top_n, gnomad, gq, ad, rankvar_filter, is_note, target, long_snp_caller, inheritance_mode, include_clinvar_report, allow_unphased_comphet, mito_report_tsv, clean_mito)
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
