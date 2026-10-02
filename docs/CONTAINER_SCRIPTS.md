# Container-resident PipeVar logic

Several PipeVar processes call project-specific scripts that are packaged in a
container rather than stored in this repository. Reviewers should treat the
container tag and command shown here as part of the biological implementation.

| Owner | Container | Container-resident entry points | Contract summary |
| --- | --- | --- | --- |
| `rankvar` | `beoungl/docker_test:rankvar_0.2.0` | `/opt/RankVar/RankVar.py`, `RankVar_nanocaller.py`, `normalize_gnomad_frequency.py` | ANNOVAR table, HPO and Phen2Gene evidence to `${sample}.rank_var.tsv`; applies GQ/AD, population-frequency and caller-specific depth controls. |
| `rankscore`, `rankscore_preannotated` | `beoungl/docker_test:rankscore_0.3.0` and `0.3.1` | `/rankscore/clinvar.sh`, `validate_preannotated_annovar_pair.py` | ANNOVAR evidence to `${sample}.rankscore_filtered.tsv` and `${sample}.clinvar.txt`; the preannotated route also validates TXT/VCF identity. |
| NGS/SNP/SV prioritization | `beoungl/docker_test:longphase_0.4.0` | `/assign_dom_or_rec*.py`, `/prio_gene_only.py`, `/clinvar_vcf_and_txt.sh`, `/rankscore_vcf_and_txt.sh`, `/rankvar_vcf_and_tsv.sh`, `/phenosv_vcf_and_tsv.sh`, `/filter_phenosv_vcf.py` | Ranked evidence plus phenotype to `${sample}.prio.vcf`, `${sample}.prio_gene.vcf`, and `${sample}.frequency_audit.tsv`. One age-capable process serves single-sample and CSV inputs; single inputs pass a blank age. |
| LongPhase prioritization support | `beoungl/docker_test:longphase_0.4.0` | the same prioritization scripts above | Adds phased SNV/SV evidence. The same optional age field is passed in single and CSV modes. |
| `mito_prio` | `beoungl/docker_test:mito_annotation_0.4.2` | `/opt/mito/bin/prioritize_mito_variants.py` | Annotated mitochondrial TSV/VCF evidence to `${sample}.mito.prioritized.tsv` using configured depth, allele-fraction, APOGEE2 and MitoTip thresholds. |

Some current tags are versioned, but none of the custom image references above
are content-addressed. Before changing one,
archive or vendor the invoked script sources, record the image digest, and run
variant-concordance tests. A container rebuild can change biological results
even when the Nextflow module is unchanged.

## Host resources and runtime profiles

`annovar` and `annovar_sv` mount `--annovar_host_path` at `/annovar`.
`phenosv` mounts `--phenosv_host_path` at `/PhenoSV/train_data`.
Each consuming task checks that its host directory exists before constructing
the mount. Other tasks do not require these directories. The ANNOVAR directory
must contain the installed scripts and `humandb`; the PhenoSV directory must
contain the resources prepared by `setup.sh`.

PhenoGPT2 mounts its configured model directories read-only and its optional
cache read-write. Mount arguments are quoted to support spaces and apostrophes.
All four existing profiles retain their executor and GPU backend selections.

`setup.sh` downloads resources, then calls `update_bind_paths.sh` and
`update_profile.sh`. The same helpers can update configuration later without
downloading resources again. Resource paths are stored as absolute paths.

The local `docker_work` source tree is useful for reviewing script contracts;
it does not prove what is installed in a published image. This cleanup keeps
the existing image references and does not rebuild containers. Image execution,
GPU availability, and complete database/model contents require runtime checks.

## Docker source audit notes

- `docker_work/html_report/Dockerfile` installs `python3-pip`, although its
  report script currently uses only the Python standard library and installs no
  packages with pip. Remove it only with a rebuilt-image smoke test.
- The reviewed PhenoSV and RankVar Dockerfiles install broad build and
  container-management packages. Their necessity cannot be established from
  the Nextflow commands alone, so package removal requires image build and
  representative runtime tests.
- Several custom images still use mutable tags. The local Dockerfiles and
  scripts therefore remain audit aids rather than proof of the published image
  contents; record image digests before changing or releasing them.
