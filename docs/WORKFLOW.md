# PipeVar_mito Workflow

This page summarizes PipeVar_mito as a paper-style workflow. The first figures
emphasize analytical stages and evidence flow rather than command-line routing
details. Single-sample and batch execution use the same biological workflow;
batch mode runs the corresponding `multi_*` modules per sample.

## Paper-Style Workflow Summary

```mermaid
flowchart TD
    samples["Patient sample data\nBAM/CRAM alignments or VCF files"]
    reference["Reference genome resources\nFASTA indexes, annotation databases, repeat catalog"]
    phenotype["Clinical phenotype input\nClinical note or HPO terms"]
    phenotypePrep["Phenotype processing\nphenotagger when notes are supplied\nPhen2gene gene ranking"]
    discovery["Variant discovery and normalization"]
    snv["SNV/indel evidence"]
    sv["Structural variant evidence"]
    cnv["Copy-number variant evidence"]
    mei["Mobile element insertion evidence"]
    repeat["Repeat expansion evidence"]
    mito["Mitochondrial variant evidence"]
    annotation["Variant annotation and evidence integration\nANNOVAR, mtDNA annotations, phenotype scores"]
    prioritization["Disease-gene prioritization\nRankVar, RankScore, PhenoSV, ngs_prio or longphase"]
    report["Final candidate reports\nprioritized VCF/TSV outputs and HTML report"]

    samples --> discovery
    reference --> discovery
    reference --> annotation
    phenotype --> phenotypePrep --> annotation
    phenotypePrep --> prioritization
    discovery --> snv
    discovery --> sv
    discovery --> cnv
    discovery --> mei
    discovery --> repeat
    discovery --> mito
    snv --> annotation
    sv --> annotation
    cnv --> annotation
    mei --> annotation
    repeat --> prioritization
    mito --> annotation
    annotation --> prioritization --> report
```

## Nuclear Variant Analysis

```mermaid
flowchart TD
    input["Nuclear input data\nshort-read, long-read, or pre-called VCF"]
    phenotype["Phenotype evidence\nHPO terms from input or phenotagger"]
    phen2gene["Phen2gene\nphenotype-ranked gene list"]
    target["Optional phenotype-targeted regions\nreduce_region_phen2gene"]

    subgraph snvLane["SNV/indel analysis"]
        snvCall["SNV/indel calling\nshort reads: DeepVariant or HaplotypeCaller\nlong reads: Clair3 or NanoCaller\nVCF input: use supplied calls"]
        snvAnno["ANNOVAR annotation"]
        rankvar["RankVar rare-disease scoring"]
        rankscore["RankScore and ClinVar evidence filtering"]
    end

    subgraph svLane["Structural variant, CNV, and mobile element analysis"]
        svCall["Breakpoint/SV calling\nshort reads: Manta\nlong reads: Sniffles\nVCF input: use supplied calls"]
        cnvCall["Copy-number variant detection\nshort reads: CNVnator\nlong reads: CNVpytor"]
        meiCall["Mobile element insertion detection\nshort reads: xTEA"]
        svMerge["SV/CNV/MEI evidence merge\nmerge_shortread_sv_callers or merge_longread_sv_callers"]
        svAnno["ANNOVAR_SV annotation"]
        survivor["SURVIVOR consolidation"]
        phenosv["PhenoSV phenotype-aware SV scoring"]
    end

    subgraph repeatLane["Repeat expansion analysis"]
        eh["ExpansionHunter for short-read BAM/CRAM"]
        ehFilter["eh_filter"]
        nanoRepeat["NanoRepeat for long-read BAM/CRAM"]
    end

    subgraph integrate["Integrated nuclear prioritization"]
        ngsPrio["ngs_prio\nshort-read combined prioritization"]
        longphase["longphase\nlong-read haplotype-aware prioritization"]
        snpOnly["snp_prio for SNP-only workflows"]
        svOnly["sv_prio for SV-only workflows"]
    end

    output["Prioritized nuclear outputs\n*.prio.vcf, *.prio_gene.vcf, evidence tables, HTML report"]

    input --> snvCall
    input --> svCall
    input --> cnvCall
    input --> meiCall
    input --> eh
    input --> nanoRepeat
    phenotype --> phen2gene
    phen2gene --> target
    phen2gene --> snvAnno
    phen2gene --> svAnno
    target -.-> snvCall
    target -.-> snvAnno
    target -.-> svAnno

    snvCall --> snvAnno --> rankvar --> ngsPrio
    snvAnno --> rankscore --> ngsPrio
    rankvar --> longphase
    rankscore --> longphase
    rankvar --> snpOnly
    rankscore --> snpOnly

    svCall --> svMerge
    cnvCall --> svMerge
    meiCall --> svMerge
    svMerge --> svAnno --> survivor --> phenosv --> ngsPrio
    phenosv --> longphase
    phenosv --> svOnly

    eh --> ehFilter --> ngsPrio
    nanoRepeat --> longphase

    ngsPrio --> output
    longphase --> output
    snpOnly --> output
    svOnly --> output
```

## Mitochondrial Variant Calling and Annotation

```mermaid
flowchart TD
    input["BAM/CRAM alignment with index"]
    contig["Mitochondrial contig selection\npreferred alias with fallback to MT, chrM, M, or chrMT"]

    subgraph shortRead["Short-read mtDNA calling"]
        subset["Extract mitochondrial reads"]
        rg["Add read group if missing"]
        revert["RevertSam to unmapped BAM"]
        fastq["SamToFastq"]
        realign["BWA MEM realignment to reference"]
        merge["MergeBamAlignment"]
        markdup["MarkDuplicates"]
        mutect["Mutect2 mitochondria mode"]
        filter["FilterMutectCalls"]
        shortNorm["bcftools normalization"]
    end

    subgraph longRead["Long-read mtDNA calling"]
        clair3["Clair3 haploid-sensitive calling on mtDNA contig"]
        longRaw["Raw Clair3 mtDNA VCF\n*.mito.clair3.raw.vcf.gz"]
        longNorm["bcftools normalization"]
        adapt["clair3_mito_adapt.py\nadapt long-read VCF to mtDNA annotation contract"]
    end

    mitoVcf["Normalized mtDNA VCF\n*.mito.vcf.gz"]

    subgraph annotate["mtDNA annotation and prioritization"]
        annotation["annotate_mito_variants.py\nMITOMAP, MitoTip, t-APOGEE, MitImpact, HmtVar status"]
        annotated["Annotated mtDNA outputs\n*.mito.annotated.tsv\n*.mito.annotated.vcf.gz"]
        prioritize["prioritize_mito_variants.py\nheteroplasmy, depth, alternate-read, APOGEE2 and MitoTip filters"]
        prioritized["Prioritized mtDNA table\n*.mito.prioritized.tsv"]
    end

    report["Mitochondrial rows available to final HTML report\nwhen nuclear full workflow is also run"]

    input --> contig
    contig --> subset
    contig --> clair3
    subset --> rg --> revert --> fastq --> realign --> merge --> markdup --> mutect --> filter --> shortNorm --> mitoVcf
    clair3 --> longRaw --> longNorm --> adapt --> mitoVcf
    mitoVcf --> annotation --> annotated --> prioritize --> prioritized --> report
```

## Reading the Figures

- These first diagrams describe analytical stages, not every run-mode branch
  from `main.nf`.
- Short-read and long-read paths differ mainly at the caller and integration
  steps: short-read full analysis uses `ngs_prio`, while long-read full analysis
  uses `longphase`.
- VCF input bypasses read-level calling and enters the relevant SNP or SV
  annotation/prioritization lane.
- CSV mode follows the same scientific workflow as single-sample mode, using
  batch-aware `multi_*` modules.
- CNV detection and mobile element insertion detection are represented as
  distinct evidence sources in the nuclear workflow. In the current
  implementation, short-read CNVs are called with CNVnator, long-read CNVs with
  CNVpytor, and short-read mobile element insertions with xTEA.
- Optional capabilities such as targeted calling and mitochondrial analysis are
  shown as analysis stages rather than central routing decisions.

## Implementation Routing Appendix

The following diagrams retain implementation-level routing details for readers
who need to map the paper-style workflow back to `main.nf`.

### Top-Level Routing

```mermaid
flowchart TD
    start["Inputs + params"]
    validate["Normalize and validate params"]
    ref["Reference/index setup"]
    pheno["Phenotype source: note or HPO"]
    inputMode{"Input mode"}
    csv["CSV batch wrapper"]
    single["Single-sample wrapper"]
    csvAge["Optional CSV age_of_onset/age"]
    dataType{"Data type"}
    vcf["VCF mode"]
    bam["BAM/CRAM mode"]
    readType{"--type"}
    short["Short-read BAM/CRAM"]
    long["Long-read BAM/CRAM: ont or pacbio"]
    vcfMode{"--mode"}
    shortMode{"--mode"}
    longMode{"--mode"}
    vcfSnp["VCF SNP re-annotation/prioritization"]
    vcfSv["VCF SV re-annotation/prioritization"]
    ngsSnp["Short-read SNP path"]
    ngsSv["Short-read SV path"]
    ngsAll["Short-read full path: SNP + SV + repeat"]
    longSnp["Long-read SNP path"]
    longSv["Long-read SV path"]
    longAll["Long-read full path: SNP + SV + repeat + longphase"]
    shortMitoGuard{"--mito yes?"}
    longMitoGuard{"--mito yes?"}
    mitoShort["Short-read mito add-on"]
    mitoLong["Long-read mito add-on"]
    reports["Prioritized outputs and HTML report"]

    start --> validate --> ref --> pheno --> inputMode
    inputMode -->|"--input_csv"| csv --> csvAge --> dataType
    inputMode -->|"single --bam or --vcf"| single --> dataType
    dataType -->|"--vcf true or --vcf file"| vcf --> vcfMode
    dataType -->|"--bam true or --bam file"| bam --> readType
    readType -->|"short"| short --> shortMitoGuard
    readType -->|"ont or pacbio"| long --> longMitoGuard
    vcfMode -->|"snp required"| vcfSnp --> reports
    vcfMode -->|"sv required"| vcfSv --> reports
    shortMitoGuard -->|"yes, snp/full only"| mitoShort
    shortMitoGuard -->|"no or not allowed"| shortMode
    longMitoGuard -->|"yes, snp/full only, not light"| mitoLong
    longMitoGuard -->|"no or not allowed"| longMode
    mitoShort -.-> ngsAll
    mitoLong -.-> longAll
    shortMode -->|"snp"| ngsSnp --> reports
    shortMode -->|"sv"| ngsSv --> reports
    shortMode -->|"omitted"| ngsAll --> reports
    longMode -->|"snp"| longSnp --> reports
    longMode -->|"sv"| longSv --> reports
    longMode -->|"omitted"| longAll --> reports
```

### Short-Read BAM/CRAM Path

```mermaid
flowchart TD
    input["Short-read BAM/CRAM + index"]
    ref["Reference FASTA + .fai"]
    phenoChoice{"Phenotype input"}
    note["Clinical note"]
    hpoFile["HPO IDs"]
    tagger["phenotagger"]
    hpo["HPO terms"]
    p2g["Phen2gene"]
    targetChoice{"--target yes?"}
    targetBed["reduce_region_phen2gene target BED"]
    snpCaller{"SNP caller"}
    deepvariant["deepvariant default"]
    haplotypecaller["haplotypecaller with --light yes"]
    annovar["ANNOVAR"]
    rankvar["RankVar"]
    rankscore["Rankscore_analysis"]
    svCaller["Manta"]
    xteaChoice{"--xtea yes?"}
    xtea["xTEA"]
    cnvnatorChoice{"--cnvnator yes?"}
    normalize["normalize_shortread_alignment"]
    cnvnator["CNVnator"]
    mergeSv["merge_shortread_sv_callers"]
    annovarSv["ANNOVAR_SV"]
    survivor["SURVIVOR"]
    phenosv["PhenoSV"]
    expansion["ExpansionHunter"]
    ehFilter["eh_filter"]
    ngsPrio["ngs_prio"]
    html{"--mito yes?"}
    report["variant_html_report"]
    reportMito["variant_html_report_with_mito"]
    mitoPrep["mito_prep_mutect2"]
    mitoMutect["mito_mutect2"]
    mitoAnno["mito_annotation"]
    mitoPrio["mito_prio"]

    input --> snpCaller
    ref --> snpCaller
    phenoChoice -->|"--note or CSV note mode"| note --> tagger --> hpo
    phenoChoice -->|"--hpo or CSV --note no"| hpoFile --> hpo
    hpo --> p2g --> targetChoice
    targetChoice -->|"yes"| targetBed --> snpCaller
    targetChoice -->|"no"| snpCaller
    snpCaller --> deepvariant --> annovar
    snpCaller --> haplotypecaller --> annovar
    p2g --> annovar
    annovar --> rankvar --> ngsPrio
    annovar --> rankscore --> ngsPrio

    input --> svCaller --> mergeSv
    input --> xteaChoice
    ref --> xtea
    xteaChoice -->|"yes, SV/full only"| xtea --> mergeSv
    input --> cnvnatorChoice
    cnvnatorChoice -->|"yes by default, SV/full only"| normalize --> cnvnator --> mergeSv
    mergeSv --> annovarSv --> survivor --> phenosv --> ngsPrio

    input --> expansion --> ehFilter --> report
    ehFilter --> reportMito
    ngsPrio --> html
    html -->|"no"| report
    html -->|"yes"| reportMito

    input -.->|"--mito yes, SNP/full only"| mitoPrep --> mitoMutect --> mitoAnno --> mitoPrio --> reportMito
```

### Long-Read BAM/CRAM Path

```mermaid
flowchart TD
    input["Long-read BAM/CRAM + index"]
    ref["Reference FASTA + .fai"]
    phenoChoice{"Phenotype input"}
    note["Clinical note"]
    hpoFile["HPO IDs"]
    tagger["phenotagger"]
    hpo["HPO terms"]
    p2g["Phen2gene"]
    targetChoice{"--target yes?"}
    targetBed["reduce_region_phen2gene target BED"]
    snpCaller{"SNP caller"}
    clair3["clair3 default"]
    nanocaller["nanocaller with --light yes"]
    annovar["ANNOVAR"]
    rankvar["RankVar"]
    rankscore["Rankscore_analysis"]
    sniffles["Sniffles"]
    cnvpytorChoice{"--cnvpytor yes?"}
    cnvpytor["CNVpytor"]
    mergeSv["merge_longread_sv_callers"]
    annovarSv["ANNOVAR_SV"]
    survivor["SURVIVOR"]
    phenosv["PhenoSV"]
    repeat["NanoRepeat"]
    longphase["longphase"]
    html{"--mito yes?"}
    report["variant_html_report"]
    reportMito["variant_html_report_with_mito"]
    mitoClair3["mito_clair3"]
    mitoClair3Post["mito_clair3_postprocess"]
    mitoAnno["mito_annotation"]
    mitoPrio["mito_prio"]

    input --> snpCaller
    ref --> snpCaller
    phenoChoice -->|"--note or CSV note mode"| note --> tagger --> hpo
    phenoChoice -->|"--hpo or CSV --note no"| hpoFile --> hpo
    hpo --> p2g --> targetChoice
    targetChoice -->|"yes"| targetBed --> snpCaller
    targetChoice -->|"no"| snpCaller
    snpCaller --> clair3 --> annovar
    snpCaller --> nanocaller --> annovar
    p2g --> annovar
    annovar --> rankvar --> longphase
    annovar --> rankscore --> longphase

    input --> sniffles
    input --> cnvpytorChoice
    annovar -.->|"SNP VCF for BAF in full mode"| cnvpytorChoice
    cnvpytorChoice -->|"yes, SV/full only; default no"| cnvpytor
    sniffles --> mergeSv
    cnvpytor --> mergeSv
    mergeSv --> annovarSv --> survivor --> phenosv --> longphase

    input --> repeat --> longphase
    longphase --> html
    html -->|"no"| report
    html -->|"yes"| reportMito

    input -.->|"--mito yes, SNP/full only, not --light yes"| mitoClair3 --> mitoClair3Post --> mitoAnno --> mitoPrio --> reportMito
```

### Mitochondrial Calling and Annotation Detail

```mermaid
flowchart TD
    enable{"--mito yes"}
    bam["BAM/CRAM + index"]
    refShort["Short-read reference bundle: FASTA + .fai + .dict + BWA indexes"]
    refLong["Long-read reference bundle: FASTA + .fai"]
    contig["Resolve mitochondrial contig: mito_contig, MT, chrM, M, chrMT"]
    readType{"Read type"}

    subgraph short["Short-read variant calling"]
        subset["Subset mitochondrial reads"]
        rg["Add read group if missing"]
        revert["GATK RevertSam"]
        fastq["GATK SamToFastq"]
        realign["BWA MEM realignment"]
        merge["GATK MergeBamAlignment"]
        markdup["GATK MarkDuplicates"]
        mutect["GATK Mutect2 mitochondria mode"]
        filter["GATK FilterMutectCalls"]
        normShort["bcftools norm split multiallelics"]
        shortVcf["*.mito.vcf.gz + .tbi"]
    end

    subgraph long["Long-read variant calling"]
        clair["run_clair3.sh haploid-sensitive on mito contig"]
        rawLong["*.mito.clair3.raw.vcf.gz"]
        normLong["bcftools norm split multiallelics"]
        adapt["clair3_mito_adapt.py"]
        longVcf["*.mito.vcf.gz + .tbi"]
    end

    subgraph annotate["Shared annotation and prioritization"]
        anno["annotate_mito_variants.py"]
        annoTsv["*.mito.annotated.tsv"]
        annoVcf["*.mito.annotated.vcf.gz + .tbi"]
        prio["prioritize_mito_variants.py"]
        filters["VAF/depth/alt-read and GUI score filters"]
        prioTsv["*.mito.prioritized.tsv"]
        final["Full-mode HTML report can include mito TSV"]
    end

    enable --> bam --> contig --> readType
    refShort --> subset
    refShort --> mutect
    refLong --> clair
    readType -->|"--type short"| subset
    readType -->|"--type ont/pacbio, not --light yes"| clair

    subset --> rg --> revert --> fastq --> realign --> merge --> markdup --> mutect --> filter --> normShort --> shortVcf
    clair --> rawLong --> normLong --> adapt --> longVcf
    shortVcf --> anno
    longVcf --> anno
    anno --> annoTsv --> prio
    anno --> annoVcf
    prio --> filters --> prioTsv --> final
```

### Guarded Branch Legend

- `--light yes`: changes the SNP caller only. Short reads use
  `haplotypecaller`; long reads use `nanocaller`. Defaults are `deepvariant`
  for short reads and `clair3` for long reads.
- `--target yes`: inserts `Phen2gene -> reduce_region_phen2gene`, then passes
  the target BED into SNP calling and downstream annotation.
- `--mito yes`: BAM/CRAM only, SNP/full only, never VCF-only or SV-only.
  Short-read mito uses Mutect2; long-read mito uses Clair3 and is rejected with
  long-read `--light yes`.
- `--xtea yes`: short-read BAM/CRAM SV/full only. xTEA VCFs merge with
  Manta and optional CNVnator before `ANNOVAR_SV`.
- `--cnvnator yes`: short-read SV/full branch, default on. The path normalizes
  the alignment before `CNVnator`.
- `--cnvpytor yes`: long-read SV/full branch, default off. In SV-only mode it
  runs read-depth only; in full mode it can use SNP/BAF support from the SNP VCF.

### Verification Notes

The diagrams intentionally omit legacy or undispatched paths such as the old
`*_light` subworkflows and `cuteSV`. Current `main.nf` dispatches unified
workflows and passes caller mode with `short_snp_caller` and `long_snp_caller`.

To check the diagrams against the current source:

```bash
rg -n "SINGLE_ALIGNMENT_|INPUT_CSV_|short_snp_caller|long_snp_caller|clean_mito|clean_xtea|clean_cnvpytor" main.nf
rg -n "workflow SINGLE_ALIGNMENT_ALL_NGS|workflow SINGLE_ALIGNMENT_ALL_LONGPHASE|workflow SINGLE_ALIGNMENT_NGS_MITO|workflow SINGLE_ALIGNMENT_LONG_MITO|include \\{" subworkflows/*/main.nf
rg -n "CNVnator|CNVpytor|xTEA|xtea|merge_shortread_sv_callers|merge_longread_sv_callers" subworkflows modules
```
