// Validate inputs once, then prepare sample-keyed records for every biological route.
// No analysis tasks run here. See workflow.md for the route map and record glossary.
workflow INPUT_PREPARATION {
    main:
    // ------------------------------------------------------------------
    // 0. NORMALIZE CLI OPTIONS (trim whitespace and normalize option case)
    // ------------------------------------------------------------------

    // Normalize option spelling once before validation and workflow selection.

    def clean_mode   = params.mode   ? params.mode.trim().toLowerCase()   : null
    def clean_type   = params.type   ? params.type.trim().toLowerCase()   : 'ont' // Default to Oxford Nanopore when omitted.
    def clean_light  = params.light == null ? 'no' : params.light.toString().trim().toLowerCase()
    def clean_genome = params.genome ? params.genome.trim().toLowerCase() : 'hg38'
    def raw_target_param = params.target ?: params.targeted
    def clean_target = raw_target_param ? raw_target_param.trim().toLowerCase() : 'no'
    def clean_note   = params.note   ? params.note.trim().toLowerCase()   : 'no'
    def clean_inheritance_mode = params.inheritance_mode ? params.inheritance_mode.trim().toLowerCase() : 'ml'
    def clean_include_clinvar_report = params.include_clinvar_report ? params.include_clinvar_report.trim().toLowerCase() : 'yes'
    def clean_allow_unphased_comphet = params.allow_unphased_comphet ? params.allow_unphased_comphet.trim().toLowerCase() : 'no'
    def clean_prioritize_sv_only = params.prioritize_sv_only ? params.prioritize_sv_only.toString().trim().toLowerCase() : 'no'
    def clean_common_sv_filter = params.common_sv_filter ? params.common_sv_filter.toString().trim().toLowerCase() : 'no'
    def clean_rankscore_softwares = params.rankscore_softwares ? params.rankscore_softwares.toString().trim() : ""
    def normalized_gpu_mode = params.GPU ? params.GPU.toString().trim().toLowerCase() : 'no'
    def clean_gpu_backend = params.gpu_backend ? params.gpu_backend.toString().trim().toLowerCase() : 'singularity'
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
    def clean_xtea = params.xtea ? params.xtea.toString().trim().toLowerCase() : 'no'
    def clean_mito = params.mito ? params.mito.toString().trim().toLowerCase() : 'no'
    def clean_annotated_snv = params.annotated_snv ? params.annotated_snv.toString().trim().toLowerCase() : 'no'
    def clean_annotated_sv = params.annotated_sv ? params.annotated_sv.toString().trim().toLowerCase() : 'no'

    def gnomadAdCeilingText = params.gnomad_af_ad?.toString()?.trim()
    def gnomadArCeilingText = params.gnomad_af_ar?.toString()?.trim()

    def parseFrequencyThreshold = { String name, String text ->
        double value
        try {
            value = Double.parseDouble(text)
        }
        catch (NumberFormatException ignored) {
            error "ERROR: --${name} must be a finite numeric value within [0,1] (received '${text}')."
        }
        if (!Double.isFinite(value) || value < 0.0d || value > 1.0d) {
            error "ERROR: --${name} must be a finite numeric value within [0,1] (received '${text}')."
        }
        value
    }
    def validatedGnomadAdCeiling = parseFrequencyThreshold('gnomad_af_ad', gnomadAdCeilingText)
    def validatedGnomadArCeiling = parseFrequencyThreshold('gnomad_af_ar', gnomadArCeilingText)
    if (validatedGnomadAdCeiling > validatedGnomadArCeiling) {
        error "ERROR: --gnomad_af_ad (${gnomadAdCeilingText}) must be less than or equal to --gnomad_af_ar (${gnomadArCeilingText})."
    }
    params.gnomad_af_ad = gnomadAdCeilingText
    params.gnomad_af_ar = gnomadArCeilingText

    def deprecatedCommonSvAfText = params.common_sv_af?.toString()?.trim()
    def commonSvAdCeilingText = params.common_sv_af_ad?.toString()?.trim()
    def commonSvArCeilingText = params.common_sv_af_ar?.toString()?.trim()
    if (deprecatedCommonSvAfText && (commonSvAdCeilingText || commonSvArCeilingText)) {
        error "ERROR: --common_sv_af is deprecated and cannot be combined with --common_sv_af_ad or --common_sv_af_ar."
    }
    if (deprecatedCommonSvAfText) {
        log.warn "--common_sv_af is deprecated; applying '${deprecatedCommonSvAfText}' to both AD and AR common-SV ceilings."
        commonSvAdCeilingText = deprecatedCommonSvAfText
        commonSvArCeilingText = deprecatedCommonSvAfText
    }
    else {
        commonSvAdCeilingText = commonSvAdCeilingText ?: '0.005'
        commonSvArCeilingText = commonSvArCeilingText ?: '0.01'
    }
    def validatedCommonSvAdCeiling = parseFrequencyThreshold('common_sv_af_ad', commonSvAdCeilingText)
    def validatedCommonSvArCeiling = parseFrequencyThreshold('common_sv_af_ar', commonSvArCeilingText)
    if (validatedCommonSvAdCeiling > validatedCommonSvArCeiling) {
        error "ERROR: --common_sv_af_ad (${commonSvAdCeilingText}) must be less than or equal to --common_sv_af_ar (${commonSvArCeilingText})."
    }
    params.common_sv_af_ad = commonSvAdCeilingText
    params.common_sv_af_ar = commonSvArCeilingText

    // ------------------------------------------------------------------
    // 1. VALIDATE MANIFEST AND OPTION VOCABULARIES
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
    def valid_gpu_backends = ['singularity', 'docker']

    if (!valid_gpu_backends.contains(clean_gpu_backend)) {
        error "ERROR: Invalid --gpu_backend '${params.gpu_backend}'. Use 'singularity' or 'docker'."
    }

    def validateStagedSif = { rawValue, String option ->
        if (rawValue == null || !rawValue.toString().trim()) {
            return null
        }
        if (clean_gpu_backend != 'singularity') {
            error "ERROR: --${option} may be used only with a Singularity profile."
        }
        def rawPath = rawValue.toString().trim()
        if (java.util.regex.Pattern.compile('[\\x00-\\x1F\\x7F]').matcher(rawPath).find()) {
            error "ERROR: --${option} contains a control character."
        }
        def image = new File(rawPath)
        if (!image.isAbsolute()) {
            error "ERROR: --${option} must be an absolute path: '${rawPath}'."
        }
        def canonical = image.canonicalFile
        if (!canonical.isFile()) {
            error "ERROR: --${option} must be a pre-existing regular file: '${rawPath}'."
        }
        if (!canonical.canRead()) {
            error "ERROR: --${option} is not readable: '${canonical}'."
        }
        canonical.absolutePath
    }

    params.phenogpt2_sif = validateStagedSif(params.phenogpt2_sif, 'phenogpt2_sif')
    params.phenotagger_sif = validateStagedSif(params.phenotagger_sif, 'phenotagger_sif')

    if (!params.input_csv) {
        def singlePrefix = params.out_prefix?.toString()?.trim()
        if (!singlePrefix || !(singlePrefix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
            error "ERROR: --out_prefix must use only letters, numbers, '.', '_' and '-', and cannot be empty."
        }
    }

    // Parse once into normalized, trimmed row maps. Validation and runtime channel
    // construction intentionally consume the same records.
    def parseCsvRecord = { String line ->
        def values = []
        def current = new StringBuilder()
        def quoted = false
        def skipNext = false
        line.toCharArray().eachWithIndex { ch, idx ->
            if (skipNext) {
                skipNext = false
            }
            else if (ch == ('"' as char)) {
                if (quoted && idx + 1 < line.size() && line.charAt(idx + 1) == ('"' as char)) {
                    current.append('"')
                    skipNext = true
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
    def manifestAnnotatedRowsWithoutPreannotatedSv = 0

    if (params.input_csv) {
        def manifestFile = file(params.input_csv)
        if (!manifestFile.exists()) {
            error "ERROR: input CSV not found: ${params.input_csv}"
        }

        def manifestLines = manifestFile.readLines()
        if (manifestLines.isEmpty()) {
            error "ERROR: input CSV is empty: ${params.input_csv}"
        }
        def manifestDataLines = manifestLines.drop(1).findAll { it.trim() }
        if (manifestDataLines.isEmpty()) {
            error "ERROR: Input CSV must contain at least one sample row."
        }

        manifestHeaderColumns = parseCsvRecord(manifestLines[0])
        if (manifestHeaderColumns && manifestHeaderColumns[0].startsWith('\uFEFF')) {
            manifestHeaderColumns[0] = manifestHeaderColumns[0].substring(1)
        }
        def blankHeaderColumns = manifestHeaderColumns
        .withIndex()
        .findAll { header, idx -> !header }
        .collect { header, idx -> idx + 1 }
        if (!blankHeaderColumns.isEmpty()) {
            error "ERROR: Input CSV contains blank header name(s) at column(s): ${blankHeaderColumns.join(', ')}"
        }
        def duplicateHeaders = manifestHeaderColumns
        .groupBy { it }
        .findAll { header, occurrences -> header && occurrences.size() > 1 }
        .keySet()
        if (!duplicateHeaders.isEmpty()) {
            error "ERROR: Input CSV contains duplicate header(s): ${duplicateHeaders.join(', ')}"
        }
        manifestUsesUnifiedSchema = manifestHeaderColumns.contains('input_kind')

        if (!manifestHeaderColumns.contains('sample')) {
            error "ERROR: Input CSV is missing required header: sample"
        }
        if (!manifestUsesUnifiedSchema) {
            def missingLegacyCsvHeaders = ['file_path', 'note_path'].findAll { !manifestHeaderColumns.contains(it) }
            if (!missingLegacyCsvHeaders.isEmpty()) {
                error "ERROR: Legacy input CSV is missing required header(s): ${missingLegacyCsvHeaders.join(', ')}"
            }
        }
        // Every route publishes sample-named files, so sample keys must be present
        // and unique before either the unified or legacy manifest is materialized.
        def seenManifestSamples = [] as Set
        manifestDataLines.eachWithIndex { line, rowIdx ->
            def values = parseCsvRecord(line)
            if (values.size() != manifestHeaderColumns.size()) {
                error "ERROR: Input CSV row ${rowIdx + 2} has ${values.size()} column(s); expected ${manifestHeaderColumns.size()}."
            }
            def sample = values[manifestHeaderColumns.indexOf('sample')]
            if (!sample) {
                error "ERROR: Input CSV row ${rowIdx + 2} is missing a non-empty sample value."
            }
            if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                error "ERROR: sample '${sample}' is not filename-safe; use only letters, numbers, '.', '_' and '-'."
            }
            if (seenManifestSamples.contains(sample)) {
                error "ERROR: Input CSV contains duplicate sample '${sample}'. Samples must be unique."
            }
            seenManifestSamples << sample
            def row = [:]
            manifestHeaderColumns.eachWithIndex { column, idx -> row[column] = values[idx] }
            manifestRowsForValidation << row
        }

        if (manifestUsesUnifiedSchema) {
            def requiredHeaders = ['sample', 'input_kind', 'phenotype_path', 'phenotype_format']
            def missingHeaders = requiredHeaders.findAll { !manifestHeaderColumns.contains(it) }
            if (!missingHeaders.isEmpty()) {
                error "ERROR: Unified input CSV is missing required header(s): ${missingHeaders.join(', ')}"
            }

            manifestRowsForValidation.eachWithIndex { row, rowIdx ->
                def sample = row.sample

                def inputKind = row.input_kind?.toLowerCase()
                def phenotypeFormat = row.phenotype_format?.toLowerCase()
                if (!valid_input_kinds.contains(inputKind)) {
                    error "ERROR: Unified input CSV sample '${sample}' has invalid input_kind '${row.input_kind}'."
                }
                if (!valid_phenotype_formats.contains(phenotypeFormat)) {
                    error "ERROR: Unified input CSV sample '${sample}' has invalid phenotype_format '${row.phenotype_format}'."
                }
                if (!(row.phenotype_path ?: '').toString().trim()) {
                    error "ERROR: Unified input CSV sample '${sample}' requires phenotype_path."
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
                    }
                    else {
                        if (alignmentIndexPath) {
                            error "ERROR: Unified input CSV sample '${sample}' provides alignment_index_path without alignment_path."
                        }
                        manifestAnnotatedRowsWithoutAlignment += 1
                    }
                    if (svVcfPath) {
                        manifestAnnotatedRowsWithPreannotatedSv += 1
                    }
                    else {
                        manifestAnnotatedRowsWithoutPreannotatedSv += 1
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

                manifestInputKinds << inputKind
                manifestPhenotypeFormats << phenotypeFormat
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
            if (manifestAnnotatedRowsWithPreannotatedSv > 0 && manifestAnnotatedRowsWithoutPreannotatedSv > 0) {
                error "ERROR: Unified annotated_snv manifests must either provide sv_vcf_path for every row or leave it blank for every row in v1."
            }
        }
    }

    // A manifest's variant class and the requested analysis must agree.
    if (manifestUsesUnifiedSchema && clean_mode) {
        def manifestVcfMode = manifestInputKinds == (['vcf_snv'] as Set) ? 'snp' :
        (manifestInputKinds == (['vcf_sv'] as Set) ? 'sv' : null)
        if (manifestVcfMode && clean_mode != manifestVcfMode) {
            error "ERROR: --mode conflicts with unified input_kind; use --mode ${manifestVcfMode} or omit --mode."
        }
    }

    // PhenoGPT2 resources are required only when a clinical note will actually be
    // processed. HPO-only inputs continue to work without a GPU or model mount.
    def singleClinicalNote = !params.input_csv && params.note != null && !(clean_note in ['yes', 'no'])
    def unifiedClinicalNote = params.input_csv && manifestUsesUnifiedSchema && manifestPhenotypeFormats.contains('clinical_note')
    def legacyCsvClinicalNote = params.input_csv && !manifestUsesUnifiedSchema && !(params.note != null && clean_note == 'no')
    def phenogpt2WillRun = clean_phenotype_extractor == 'phenogpt2' && (singleClinicalNote || unifiedClinicalNote || legacyCsvClinicalNote)

    def parsePositiveInteger = { value, String option ->
        try {
            def text = value?.toString()?.trim()
            if (!text || !(text ==~ /[0-9]+/)) {
                error "ERROR: --${option} must be a positive integer; received '${value}'."
            }
            def parsed = Integer.parseInt(text)
            if (parsed <= 0) {
                error "ERROR: --${option} must be a positive integer; received '${value}'."
            }
            parsed
        }
        catch (Exception ignored) {
            error "ERROR: --${option} must be a positive integer; received '${value}'."
        }
    }

    def parseNonnegativeInteger = { value, String option ->
        try {
            def text = value?.toString()?.trim()
            if (!text || !(text ==~ /[0-9]+/)) {
                error "ERROR: --${option} must be a nonnegative integer; received '${value}'."
            }
            Integer.parseInt(text)
        }
        catch (Exception ignored) {
            error "ERROR: --${option} must be a nonnegative integer; received '${value}'."
        }
    }

    params.gpu_cpus = parsePositiveInteger(params.gpu_cpus, 'gpu_cpus')
    if (params.deepvariant_max_forks != null && params.deepvariant_max_forks.toString().trim()) {
        params.deepvariant_max_forks = parsePositiveInteger(
            params.deepvariant_max_forks,
            'deepvariant_max_forks'
        )
    }

    if (phenogpt2WillRun) {
        if (normalized_gpu_mode != 'yes') {
            error "ERROR: PhenoGPT2 will process clinical notes and requires --GPU yes."
        }
        if (!(clean_phenogpt2_negation in ['yes', 'no'])) {
            error "ERROR: --phenogpt2_negation must be yes or no; received '${params.phenogpt2_negation}'."
        }

        def phenogpt2BatchSize = parsePositiveInteger(params.phenogpt2_batch_size, 'phenogpt2_batch_size')
        def phenogpt2ChunkBatchSize = parsePositiveInteger(params.phenogpt2_chunk_batch_size, 'phenogpt2_chunk_batch_size')
        def phenogpt2MaxForks = parsePositiveInteger(params.phenogpt2_max_forks, 'phenogpt2_max_forks')
        if (phenogpt2BatchSize != 1 || phenogpt2ChunkBatchSize != 1) {
            error "ERROR: PhenoGPT2 currently requires --phenogpt2_batch_size 1 and --phenogpt2_chunk_batch_size 1."
        }
        if (phenogpt2MaxForks != 1) {
            error "ERROR: PhenoGPT2 currently requires --phenogpt2_max_forks 1."
        }
        def phenogpt2Wc = parseNonnegativeInteger(params.phenogpt2_wc, 'phenogpt2_wc')
        if (phenogpt2Wc > 0) {
            error "ERROR: --phenogpt2_wc > 0 is not supported until the BERT filtering model is provisioned."
        }

        def validateHostDirectory = { rawValue, String option, boolean writable ->
            if (rawValue == null || !rawValue.toString().trim()) {
                error "ERROR: --${option} is required when PhenoGPT2 processes clinical notes."
            }
            def rawPath = rawValue.toString()
            if (java.util.regex.Pattern.compile('[,:\\x00-\\x1F\\x7F]').matcher(rawPath).find()) {
                error "ERROR: --${option} contains a bind delimiter or control character: '${rawPath}'."
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

        def mainCheckpointRoot = validateHostDirectory(
            params.phenogpt2_model_host_path,
            'phenogpt2_model_host_path',
            false
        )
        def negationCheckpointRoot = null
        def embeddingCheckpointRoot = null
        if (clean_phenogpt2_negation == 'yes') {
            negationCheckpointRoot = validateHostDirectory(
                params.phenogpt2_negation_model_host_path,
                'phenogpt2_negation_model_host_path',
                false
            )
            embeddingCheckpointRoot = validateHostDirectory(
                params.phenogpt2_embedding_model_host_path,
                'phenogpt2_embedding_model_host_path',
                false
            )
        }
        def cacheRoot = params.phenogpt2_cache_host_path ?
            validateHostDirectory(params.phenogpt2_cache_host_path, 'phenogpt2_cache_host_path', true) :
            null

        params.phenogpt2_model_host_path = mainCheckpointRoot.absolutePath
        params.phenogpt2_negation_model_host_path = negationCheckpointRoot?.absolutePath
        params.phenogpt2_embedding_model_host_path = embeddingCheckpointRoot?.absolutePath
        params.phenogpt2_cache_host_path = cacheRoot?.absolutePath
    }

    // Validate normalized option vocabularies before selecting a route.
    if (clean_mode && !valid_modes.contains(clean_mode)) {
        error "ERROR: Invalid --mode '${params.mode}'. Use 'snp' or 'sv'."
    }
    if (!valid_types.contains(clean_type)) {
        error "ERROR: Invalid --type '${params.type}'. Use 'ont', 'pacbio', or 'short'."
    }
    if (!valid_yes_no.contains(clean_light)) {
        error "ERROR: Invalid --light '${params.light}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(normalized_gpu_mode)) {
        error "ERROR: Invalid --GPU '${params.GPU}'. Use 'yes' or 'no'."
    }
    if (!valid_phenotype_extractors.contains(clean_phenotype_extractor)) {
        error "ERROR: Invalid --phenotype_extractor '${params.phenotype_extractor}'. Use 'phenotagger' or 'phenogpt2'."
    }
    if (!valid_yes_no.contains(clean_mito)) {
        error "ERROR: Invalid --mito '${params.mito}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_annotated_snv)) {
        error "ERROR: Invalid --annotated_snv '${params.annotated_snv}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_annotated_sv)) {
        error "ERROR: Invalid --annotated_sv '${params.annotated_sv}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_cnvnator)) {
        error "ERROR: Invalid --cnvnator '${params.cnvnator}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_xtea)) {
        error "ERROR: Invalid --xtea '${params.xtea}'. Use 'yes' or 'no'."
    }

    if (clean_xtea == 'yes') {
        if (params.vcf) {
            error "ERROR: --xtea yes requires BAM/CRAM input."
        }
        if (clean_type != 'short') {
            error "ERROR: --xtea yes currently supports only --type short."
        }
        if (clean_mode == 'snp') {
            error "ERROR: --xtea yes requires --mode sv or omitted --mode."
        }
    }

    if (clean_mito == 'yes') {
        def unifiedVcfInput = manifestUsesUnifiedSchema && manifestInputKinds.any {
            it in ['vcf_snv', 'vcf_sv']
        }
        def singleRawVcfInput = params.vcf && clean_annotated_snv != 'yes'
        def singleAnnotatedWithoutAlignment = !params.input_csv &&
            clean_annotated_snv == 'yes' && !params.bam
        if (singleRawVcfInput || unifiedVcfInput || singleAnnotatedWithoutAlignment) {
            error "ERROR: --mito yes requires BAM/CRAM input."
        }
        if (clean_mode == 'sv') {
            error "ERROR: --mito yes requires --mode snp or omitted --mode."
        }
        if (clean_type != 'short' && clean_light == 'yes') {
            error "ERROR: Long-read --mito yes requires Clair3; disable --light or use --type short."
        }
    }

    if (!(params.cnvnator_bin_size.toString() ==~ /[1-9][0-9]*/)) {
        error "ERROR: Invalid --cnvnator_bin_size '${params.cnvnator_bin_size}'. Provide a positive integer."
    }

    params.mode = clean_mode
    params.type = clean_type
    params.light = clean_light
    params.genome = clean_genome
    params.cnvnator = clean_cnvnator
    params.xtea = clean_xtea
    params.mito = clean_mito
    params.annotated_snv = clean_annotated_snv
    params.annotated_sv = clean_annotated_sv

    if (!valid_genomes.contains(clean_genome)) {
        error "ERROR: Unsupported --genome '${params.genome}'. Use 'hg38' or 'grch38'."
    }
    if (!valid_inheritance_modes.contains(clean_inheritance_mode)) {
        error "ERROR: Invalid --inheritance_mode '${params.inheritance_mode}'. Use 'ml', 'omim', or 'gnomad'."
    }
    if (!valid_yes_no.contains(clean_include_clinvar_report)) {
        error "ERROR: Invalid --include_clinvar_report '${params.include_clinvar_report}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_allow_unphased_comphet)) {
        error "ERROR: Invalid --allow_unphased_comphet '${params.allow_unphased_comphet}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_prioritize_sv_only)) {
        error "ERROR: Invalid --prioritize_sv_only '${params.prioritize_sv_only}'. Use 'yes' or 'no'."
    }
    if (!valid_yes_no.contains(clean_common_sv_filter)) {
        error "ERROR: Invalid --common_sv_filter '${params.common_sv_filter}'. Use 'yes' or 'no'."
    }
    params.common_sv_filter = clean_common_sv_filter

    if (!valid_yes_no.contains(clean_target)) {
        error "ERROR: Invalid --target/--targeted '${raw_target_param}'. Use 'yes' or 'no'."
    }

    if (!(params.phenosv_score.toString() ==~ /([0-9]+([.][0-9]+)?|[.][0-9]+)/)) {
        error "ERROR: Invalid --phenosv_score '${params.phenosv_score}'. Provide a numeric threshold."
    }

    def nanocallerDpText = params.nanocaller_dp?.toString()?.trim()
    if (!(nanocallerDpText ==~ /([0-9]+([.][0-9]+)?|[.][0-9]+)/)) {
        error "ERROR: Invalid --nanocaller_dp '${params.nanocaller_dp}'. Provide a finite, non-negative numeric threshold, for example 20."
    }

    parseFrequencyThreshold('common_sv_reciprocal_overlap', params.common_sv_reciprocal_overlap.toString())

    if (!(params.common_sv_distance.toString() ==~ /[0-9]+/)) {
        error "ERROR: Invalid --common_sv_distance '${params.common_sv_distance}'. Provide an integer distance, for example 1000."
    }

    if (!(params.common_sv_ins_distance.toString() ==~ /[0-9]+/)) {
        error "ERROR: Invalid --common_sv_ins_distance '${params.common_sv_ins_distance}'. Provide an integer distance, for example 500."
    }

    parseFrequencyThreshold('common_sv_ins_identity', params.common_sv_ins_identity.toString())

    if (clean_rankscore_softwares && clean_rankscore_softwares.split(",").every { it.trim().isEmpty() }) {
        error "ERROR: Invalid --rankscore_softwares '${params.rankscore_softwares}'. Provide a comma-separated list."
    }

    // ------------------------------------------------------------------
    // VALIDATE INPUT PRESENCE AND REFERENCE DEPENDENCIES
    // ------------------------------------------------------------------

    if (!params.input_csv && !params.bam && !params.vcf) {
        error "ERROR: No input data specified; provide --input_csv, --bam, or --vcf."
    }

    if (params.input_csv && !manifestUsesUnifiedSchema && !params.bam && !params.vcf) {
        error "ERROR: Legacy --input_csv requires either --bam true or --vcf true."
    }

    if (params.input_csv && !manifestUsesUnifiedSchema && params.bam && params.vcf) {
        error "ERROR: Legacy input CSV is ambiguous; choose exactly one of --bam true or --vcf true."
    }

    if (params.input_csv && manifestUsesUnifiedSchema && (params.bam || params.vcf)) {
        error "ERROR: Unified input CSV cannot be combined with --bam or --vcf; each row's input_kind selects the path."
    }

    if (!params.input_csv && params.bam && params.vcf && clean_annotated_snv != 'yes') {
        error "ERROR: Single-file input is ambiguous; provide either --bam or --vcf unless --annotated_snv yes combines the prepared VCF with an alignment."
    }

    if (!params.input_csv && params.note && params.hpo) {
        error "ERROR: Single-file phenotype input is ambiguous; provide exactly one of --note or --hpo."
    }

    def singleRawVcfMode = !params.input_csv && params.vcf && clean_annotated_snv != 'yes'
    def legacyCsvVcfMode = params.input_csv && !manifestUsesUnifiedSchema && params.vcf
    if ((singleRawVcfMode || legacyCsvVcfMode) && !clean_mode) {
        error "ERROR: Raw VCF input requires --mode snp or --mode sv."
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
        if (clean_mito == 'yes' && manifestAnnotatedRowsWithoutAlignment > 0) {
            error "ERROR: --mito yes requires alignment_path for every annotated_snv manifest row."
        }
    }

    // Reference bundles are required for aligned reads and documented raw-VCF routes.
    def unifiedManifestNeedsReference = manifestUsesUnifiedSchema && manifestRowsForValidation.any { row ->
        (row.input_kind ?: '').toString().trim().toLowerCase() != 'annotated_snv'
    }
    def annotatedManifestNeedsReference = manifestUsesUnifiedSchema && manifestRowsForValidation.any { row ->
        (row.input_kind ?: '').toString().trim().toLowerCase() == 'annotated_snv' && row.containsKey('alignment_path') && row.alignment_path
    }
    def unifiedManifestNeedsExpansionHunter = manifestUsesUnifiedSchema && (
        manifestInputKinds.any { it in ['bam_ngs', 'cram_ngs'] } || annotatedManifestNeedsReference
    )
    def rawVcfNeedsReference = params.vcf && clean_annotated_snv != 'yes'
    def needsReference = params.bam || rawVcfNeedsReference || unifiedManifestNeedsReference || annotatedManifestNeedsReference
    if (needsReference && !params.ref_fa) {
        error """
        ERROR: Reference FASTA missing.
        BAM/CRAM and raw VCF inputs require a reference genome.
        Please specify: --ref_fa <path/to/reference.fasta>
        """
    }

    // A supplied reference must have its FASTA index; mitochondrial GATK also needs the aligner sidecars.
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
                ERROR: Reference dictionary (.dict) not found
                Mitochondrial analysis requires a sequence dictionary for GATK.
                Expected at: ${ref_dict}
                Create it with: gatk CreateSequenceDictionary -R ${params.ref_fa}
                """
            }

            def bwa_suffixes = ['amb', 'ann', 'bwt', 'pac', 'sa']
            def missing_bwa_indexes = bwa_suffixes.findAll { suffix ->
                !file("${params.ref_fa}.${suffix}").exists()
            }
            if (!missing_bwa_indexes.isEmpty()) {
                def expected_files = missing_bwa_indexes.collect { suffix -> "${params.ref_fa}.${suffix}" }.join('\n            ')
                error """
                ERROR: BWA reference sidecar files not found
                Mitochondrial analysis requires a BWA-indexed reference bundle.
                Missing:
                ${expected_files}

                Create the sidecars with:
                bwa index ${params.ref_fa}
                """
            }
        }
    }

    def is_expansionhunter_mode = clean_type == 'short' && (
    params.bam || unifiedManifestNeedsExpansionHunter
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
            or ensure the default PipeVar_simplified/data catalog exists for --genome ${clean_genome}.
            """
        }
    }

    // Build phenotype, alignment, and annotation records.

    reference_bundle = Channel.value([])
    clinical_metadata_by_sample = null
    input_bam = null
    def bam = null
    def vcf = null
    def single_out_prefix = null
    def note = null
    mitochondrial_reference_bundle = null
    gatk_reference_bundle = null
    eh_variant_catalog = null
    def annovar_txt = null
    def annovar_sv_vcf = null
    input_annotated_snv = null
    input_annotated_snv_sv = null
    input_annotated_ngs = null
    input_annotated_called_ngs = null
    input_vcf = null
    def csv_manifest_mode = null
    def csv_manifest_is_note = null
    def phenotypeFileForRow = { row, columnName ->
        def rawPath = (row[columnName] ?: '').toString().trim()
        if (!rawPath) {
            error "ERROR: sample '${row.sample}' requires a non-empty ${columnName}."
        }
        file(rawPath, checkIfExists: true)
    }
    def resolveAlignmentIndex = { alignmentFile, explicitIndexPath, sourceLabel ->
        def explicitPath = explicitIndexPath?.toString()?.trim()
        def candidatePaths = []
        if (explicitPath) {
            candidatePaths << explicitPath
        }
        if (alignmentFile.name.endsWith('.cram')) {
            candidatePaths << "${alignmentFile}.crai"
            candidatePaths << alignmentFile.toString().replaceFirst(/\.cram$/, '.crai')
        }
        else {
            candidatePaths << alignmentFile.toString().replaceFirst(/\.bam$/, '.bai')
            candidatePaths << "${alignmentFile}.bai"
        }

        def candidate = candidatePaths.unique()
            .collect { candidatePath -> file(candidatePath) }
            .find { candidateFile -> candidateFile.exists() }
        if (candidate != null) {
            return candidate
        }
        error "Index file not found for ${alignmentFile} (${sourceLabel}). Looked for: ${candidatePaths.unique().join(', ')}"
    }
    if ( params.input_csv ) {
        clinical_metadata_by_sample = Channel.fromList(manifestRowsForValidation)
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
                if (manifestAnnotatedRowsWithPreannotatedSv > 0 && manifestAnnotatedRowsWithAlignment > 0) {
                    input_annotated_ngs = Channel.fromList(manifestRowsForValidation)
                    .map { row ->
                        def bam_file = file(row.alignment_path, checkIfExists: true)
                        def bai = resolveAlignmentIndex(bam_file, row.alignment_index_path, "unified CSV sample '${row.sample}'")
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

                    csv_manifest_mode = 'annotated_all_ngs'
                }
                else if (manifestAnnotatedRowsWithPreannotatedSv > 0) {
                    input_annotated_snv_sv = Channel.fromList(manifestRowsForValidation)
                    .map { row ->
                        return tuple(
                        row.sample,
                        file(row.snv_txt_path, checkIfExists: true),
                        file(row.snv_vcf_path, checkIfExists: true),
                        file(row.sv_vcf_path, checkIfExists: true),
                        phenotypeFileForRow(row, 'phenotype_path'),
                        (row.phenotype_format ?: '').toString().trim().toLowerCase()
                        )
                    }

                    csv_manifest_mode = 'annotated_snv_sv'
                }
                else if (manifestAnnotatedRowsWithAlignment > 0) {
                    input_annotated_called_ngs = Channel.fromList(manifestRowsForValidation)
                    .map { row ->
                        def bam_file = file(row.alignment_path, checkIfExists: true)
                        def bai = resolveAlignmentIndex(bam_file, row.alignment_index_path, "unified CSV sample '${row.sample}'")
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

                    csv_manifest_mode = 'annotated_snv_called_sv_ngs'
                }
                else {
                    input_annotated_snv = Channel.fromList(manifestRowsForValidation)
                    .map { row ->
                        return tuple(
                        row.sample,
                        file(row.snv_txt_path, checkIfExists: true),
                        file(row.snv_vcf_path, checkIfExists: true),
                        phenotypeFileForRow(row, 'phenotype_path'),
                        (row.phenotype_format ?: '').toString().trim().toLowerCase()
                        )
                    }

                    csv_manifest_mode = 'annotated_snv'
                }
            }
            else if ( manifestInputKinds == (['vcf_snv'] as Set) || manifestInputKinds == (['vcf_sv'] as Set) ) {
                input_vcf = Channel.fromList(manifestRowsForValidation)
                .map { row ->
                    return tuple(
                    row.sample,
                    file(row.vcf_path, checkIfExists: true),
                    phenotypeFileForRow(row, 'phenotype_path'),
                    )
                }

                csv_manifest_mode = manifestInputKinds.contains('vcf_sv') ? 'sv' : 'snp'
                if (manifestPhenotypeFormats.size() != 1) {
                    error "ERROR: Mixed phenotype_format values are not yet supported for non-annotated CSV modes in v1."
                }
                csv_manifest_is_note = manifestPhenotypeFormats.contains('clinical_note') ? "yes" : "no"
            }
            else {
                input_bam = Channel.fromList(manifestRowsForValidation)
                .map { row ->
                    def bam_file = file(row.alignment_path, checkIfExists: true)
                    def bai = resolveAlignmentIndex(bam_file, row.alignment_index_path, "unified CSV sample '${row.sample}'")
                    return tuple(
                    row.sample,
                    bam_file,
                    bai,
                    phenotypeFileForRow(row, 'phenotype_path'),
                    )
                }

                if (manifestPhenotypeFormats.size() != 1) {
                    error "ERROR: Mixed phenotype_format values are not yet supported for non-annotated CSV alignment modes in v1."
                }
                csv_manifest_is_note = manifestPhenotypeFormats.contains('clinical_note') ? "yes" : "no"
            }
        }
        else if ( params.vcf ) {
            input_vcf = Channel.fromList(manifestRowsForValidation)
            .map { row ->
                def vcf_file = file(row.file_path, checkIfExists: true)
                def note_file = phenotypeFileForRow(row, 'note_path')

                // Preserve the sample key beside its VCF and phenotype source.
                return tuple(
                row.sample,
                vcf_file,
                note_file,
                )
            }

        }
        else if ( params.bam ) {
            input_bam = Channel.fromList(manifestRowsForValidation)
            .map { row ->
                def bam_file = file(row.file_path, checkIfExists: true)

                def bai = resolveAlignmentIndex(bam_file, null, "legacy CSV sample '${row.sample}'")

                def note_file = phenotypeFileForRow(row, 'note_path')
                def out_prefix = row.sample

                // Alignment record: sample, aligned reads, index, phenotype source.
                return tuple(
                out_prefix,
                bam_file,
                bai,
                note_file,
                )
            }

        }
    }
    else if ( params.bam != null ) {
        bam = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { file ->

            def index = resolveAlignmentIndex(file, null, "single-sample input")

            return [ file, index ]
        }
        single_out_prefix=Channel.value(params.out_prefix)
    }
    else if ( params.vcf != null ) {
        vcf=Channel.value(file(params.vcf, checkIfExists: true))
        single_out_prefix=Channel.value(params.out_prefix)
    }
    if (clean_annotated_snv == 'yes' && !params.input_csv) {
        annovar_txt = Channel.value(file(params.annovar_txt, checkIfExists: true))
        vcf = Channel.value(file(params.vcf, checkIfExists: true))
        if (clean_annotated_sv == 'yes') {
            annovar_sv_vcf = Channel.value(file(params.annovar_sv_vcf, checkIfExists: true))
        }
    }
    if (params.ref_fa != null) {
        reference_bundle = Channel
        .fromPath(params.ref_fa, checkIfExists: true)
        .map { fa_file ->
            def fai_file = file("${fa_file}.fai")
            return [ fa_file, fai_file ]
        }
        .first()
    }
    if (params.ref_fa != null && clean_type == 'short') {
        gatk_reference_bundle = Channel
        .fromPath(params.ref_fa, checkIfExists: true)
        .map { fa_file ->
            def fai_file = file("${fa_file}.fai")
            def dict_file = file("${fa_file.parent}/${fa_file.baseName}.dict")
            return [ fa_file, fai_file, dict_file ]
        }
        .first()
    }
    if (params.ref_fa != null && clean_mito == 'yes' && clean_type == 'short') {
        mitochondrial_reference_bundle = Channel
        .fromPath(params.ref_fa, checkIfExists: true)
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
    // Clinical notes require HPO extraction; supplied HPO files enter gene ranking directly.
    def phenotype_is_clinical_note = "no"
    if ( params.input_csv ) {
        // CSV mode default: note_path is treated as clinical notes unless user explicitly sets --note no.
        phenotype_is_clinical_note = csv_manifest_is_note != null ? csv_manifest_is_note : ((params.note != null && clean_note == 'no') ? "no" : "yes")
    }
    else {
        // Single-file mode: `note` must be a file path; `hpo` is used directly as HPO input.
        if ( params.note != null && clean_note != 'yes' && clean_note != 'no' ) {
            note=Channel.value(file(params.note, checkIfExists: true))
            phenotype_is_clinical_note = "yes"
        }
        else if ( params.hpo != null ) {
            note=Channel.value(file(params.hpo, checkIfExists: true))
            phenotype_is_clinical_note = "no"
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
    mito_contig = Channel.value(params.mito_contig)
    // Single-sample scripts use OMIM for ml/omim and LOEUF for gnomad.
    def inheritance_mode_script = (clean_inheritance_mode == 'gnomad') ? 'LOEUF' : 'OMIM'
    inheritance_mode = Channel.value(inheritance_mode_script)
    include_clinvar_report = Channel.value(clean_include_clinvar_report)
    allow_unphased_comphet = Channel.value(clean_allow_unphased_comphet)
    // Use the broader AR ceiling upstream; final ranking applies AD/AR-specific ceilings.
    gnomad_af_ceiling = Channel.value(gnomadArCeilingText)
    rankscore_filter = Channel.value(params.rankscore)
    rankscore_softwares = Channel.value(clean_rankscore_softwares)
    rankvar_filter = Channel.value(params.rankvar)
    phen2gene_top_n = Channel.value(params.phen2gene_filter)
    minimum_genotype_quality = Channel.value(params.gq)
    minimum_allele_depth = Channel.value(params.ad)
    def short_read_small_variant_caller = clean_light == 'yes' ? 'haplotypecaller' : 'deepvariant'
    def long_read_small_variant_caller = clean_light == 'yes' ? 'nanocaller' : 'clair3'
    def phenotype_targeting_enabled = clean_target
    def batch_input = params.input_csv ? true : false
    def effective_mode = clean_mode ?: csv_manifest_mode
    if (!params.input_csv) {
        clinical_metadata_by_sample = single_out_prefix.map { sample -> tuple(sample, "") }
        input_bam = bam != null
        ? bam.combine(single_out_prefix).combine(note).map { bam_file, bai_file, sample, phenotype_file ->
            tuple(sample, bam_file, bai_file, phenotype_file)
        }
        : null
        input_vcf = vcf != null
        ? single_out_prefix.combine(vcf).combine(note).map { sample, vcf_file, phenotype_file ->
            tuple(sample, vcf_file, phenotype_file)
        }
        : null
        def phenotype_format = phenotype_is_clinical_note == "yes" ? "clinical_note" : "hpo"

        if (clean_annotated_snv == "yes") {
            if (params.bam != null) {
                if (clean_annotated_sv == "yes") {
                    input_annotated_ngs = single_out_prefix.combine(annovar_txt).combine(vcf).combine(annovar_sv_vcf).combine(bam).combine(note).map {
                        sample, txt_file, snv_vcf, sv_vcf, bam_file, bai_file, phenotype_file ->
                        tuple(sample, txt_file, snv_vcf, sv_vcf, bam_file, bai_file, phenotype_file, phenotype_format)
                    }
                }
                else {
                    input_annotated_called_ngs = single_out_prefix.combine(annovar_txt).combine(vcf).combine(bam).combine(note).map {
                        sample, txt_file, snv_vcf, bam_file, bai_file, phenotype_file ->
                        tuple(sample, txt_file, snv_vcf, bam_file, bai_file, phenotype_file, phenotype_format)
                    }
                }
            }
            else if (clean_annotated_sv == "yes") {
                input_annotated_snv_sv = single_out_prefix.combine(annovar_txt).combine(vcf).combine(annovar_sv_vcf).combine(note).map {
                    sample, txt_file, snv_vcf, sv_vcf, phenotype_file ->
                    tuple(sample, txt_file, snv_vcf, sv_vcf, phenotype_file, phenotype_format)
                }
            }
            else {
                input_annotated_snv = single_out_prefix.combine(annovar_txt).combine(vcf).combine(note).map {
                    sample, txt_file, snv_vcf, phenotype_file ->
                    tuple(sample, txt_file, snv_vcf, phenotype_file, phenotype_format)
                }
            }
        }
    }

    mito_input_bam = input_bam
    // Route decisions use scalar settings; records continue through named channels.
    def prepared_route = 'alignment'
    if (input_vcf != null) {
        prepared_route = 'vcf'
    }
    if (input_annotated_snv != null) {
        prepared_route = 'annotated_snv'
    }
    if (input_annotated_called_ngs != null) {
        prepared_route = 'annotated_called_ngs'
    }
    if (input_annotated_snv_sv != null) {
        prepared_route = 'annotated_snv_sv'
    }
    if (input_annotated_ngs != null) {
        prepared_route = 'annotated_ngs'
    }
    route_settings = [
    route: prepared_route,
    batch_input: batch_input,
    effective_mode: effective_mode,
    sequencing_type: clean_type,
    mitochondrial_enabled: clean_mito,
    phenotype_is_clinical_note: phenotype_is_clinical_note,
    phenotype_targeting_enabled: phenotype_targeting_enabled,
    short_read_small_variant_caller: short_read_small_variant_caller,
    long_read_small_variant_caller: long_read_small_variant_caller,
    ]
    clinical_metadata_by_sample = clinical_metadata_by_sample != null ? clinical_metadata_by_sample : Channel.empty()
    input_bam = input_bam != null ? input_bam : Channel.empty()
    input_vcf = input_vcf != null ? input_vcf : Channel.empty()
    input_annotated_snv = input_annotated_snv != null ? input_annotated_snv : Channel.empty()
    input_annotated_snv_sv = input_annotated_snv_sv != null ? input_annotated_snv_sv : Channel.empty()
    input_annotated_ngs = input_annotated_ngs != null ? input_annotated_ngs : Channel.empty()
    input_annotated_called_ngs = input_annotated_called_ngs != null ? input_annotated_called_ngs : Channel.empty()
    mito_input_bam = mito_input_bam != null ? mito_input_bam : Channel.empty()
    gatk_reference_bundle = gatk_reference_bundle != null ? gatk_reference_bundle : Channel.value([])
    mitochondrial_reference_bundle = mitochondrial_reference_bundle != null ? mitochondrial_reference_bundle : Channel.value([])
    eh_variant_catalog = eh_variant_catalog != null ? eh_variant_catalog : Channel.value([])

    emit:
    route_settings
    clinical_metadata_by_sample
    input_bam
    input_vcf
    input_annotated_snv
    input_annotated_snv_sv
    input_annotated_ngs
    input_annotated_called_ngs
    mito_input_bam
    reference_bundle
    gatk_reference_bundle
    mitochondrial_reference_bundle
    eh_variant_catalog
    mito_contig
    inheritance_mode
    include_clinvar_report
    allow_unphased_comphet
    gnomad_af_ceiling
    rankscore_filter
    rankscore_softwares
    rankvar_filter
    phen2gene_top_n
    minimum_genotype_quality
    minimum_allele_depth
}
