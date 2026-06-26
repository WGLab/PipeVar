# PipeVar Nuclear Workflow Figures

This page provides focused figure source for the completed nuclear rare-disease
analysis workflow. The diagrams emphasize analysis stages and tool handoffs,
without command-line routing details.

Rendered SVG assets:

- `docs/pipevar_nuclear_workflow.svg`
- `docs/pipevar_prioritization_workflow.svg`

## Figure 1. Nuclear Rare-Disease Workflow

```mermaid
flowchart LR
    input["Input data\nBAM/CRAM or VCF"]
    reference["Reference resources\nFASTA, annotation data, repeat catalog"]
    phenotype["Phenotype input\nClinical note or HPO terms"]

    subgraph phenotypeLane["Phenotype extraction"]
        phenotagger["PhenoTagger"]
        phen2gene["Phen2Gene"]
    end

    subgraph snvLane["SNV/indel calling"]
        snvCall["DeepVariant or HaplotypeCaller\nClair3 or NanoCaller\nSupplied VCF"]
        snvAnno["ANNOVAR"]
        snvScore["RankVar\nRankScore\nClinVar"]
    end

    subgraph svLane["SV calling"]
        svCall["Manta\nSniffles\nSupplied VCF"]
        svAnno["ANNOVAR_SV"]
        svScore["SURVIVOR\nPhenoSV"]
    end

    subgraph repeatLane["Repeat expansion analysis"]
        repeatCall["ExpansionHunter\neh_filter\nNanoRepeat"]
    end

    subgraph prioLane["Variant prioritization"]
        fullPrio["ngs_prio\nlongphase"]
        modePrio["snp_prio\nsv_prio"]
    end

    output["Final outputs\nprioritized VCF\ngene-level prioritized VCF\nrepeat table\nHTML report"]

    phenotype --> phenotagger --> phen2gene
    phenotype --> phen2gene
    input --> snvCall
    input --> svCall
    input --> repeatCall
    reference --> snvCall
    reference --> svCall
    reference --> repeatCall
    reference --> snvAnno
    reference --> svAnno
    phen2gene --> snvAnno
    phen2gene --> svAnno
    phen2gene --> fullPrio
    phen2gene --> modePrio
    snvCall --> snvAnno --> snvScore --> fullPrio
    snvScore --> modePrio
    svCall --> svAnno --> svScore --> fullPrio
    svScore --> modePrio
    repeatCall --> fullPrio
    fullPrio --> output
    modePrio --> output
```

## Figure 2. Variant Prioritization Detail

```mermaid
flowchart TD
    hpo["HPO terms"]
    inheritance["--inheritance_mode\nml | omim | gnomad"]

    subgraph snvEvidence["SNV/indel evidence"]
        annovar["ANNOVAR VCF/TXT"]
        clinvar["ClinVar"]
        rankscore["RankScore"]
        rankvar["RankVar"]
    end

    subgraph svEvidence["SV evidence"]
        annovarSv["ANNOVAR_SV VCF"]
        survivor["SURVIVOR"]
        phenosv["PhenoSV"]
    end

    subgraph assign["Inheritance-aware assignment"]
        snvAssign["assign_dom_or_rec_snp_only"]
        svAssign["assign_dom_or_rec_sv_only"]
        combinedAssign["assign_dom_or_rec"]
    end

    subgraph workflows["Prioritization workflows"]
        snpPrio["snp_prio\nSNV-only"]
        svPrio["sv_prio\nSV-only"]
        ngsPrio["ngs_prio\nshort-read full"]
        longphase["longphase\nlong-read full"]
    end

    reports["Final reports\nvariant-level VCF\ngene-level VCF"]

    annovar --> clinvar --> snvAssign
    annovar --> rankscore --> snvAssign
    annovar --> rankvar --> snvAssign
    annovarSv --> survivor --> phenosv --> svAssign
    clinvar --> combinedAssign
    rankscore --> combinedAssign
    rankvar --> combinedAssign
    phenosv --> combinedAssign
    hpo --> snvAssign
    hpo --> svAssign
    hpo --> combinedAssign
    inheritance --> snvAssign
    inheritance --> svAssign
    inheritance --> combinedAssign
    snvAssign --> snpPrio
    svAssign --> svPrio
    combinedAssign --> ngsPrio
    combinedAssign --> longphase
    snpPrio --> reports
    svPrio --> reports
    ngsPrio --> reports
    longphase --> reports
```
