# Runtime, setup, and resume

Use this reference for environment selection, resource setup, GPU/PhenoGPT2 runs, concurrency, and restart decisions. [INSTALLATION.md](../../../../docs/INSTALLATION.md) remains authoritative for detailed setup and resource guidance.

## Select an explicit profile

| Profile | Executor | Container runtime | Intended host |
| --- | --- | --- | --- |
| `local_docker` | Local | Docker | Docker-capable Linux workstation/server |
| `local_singularity` | Local | Singularity | Linux workstation/server |
| `slurm_singularity` | SLURM | Singularity | Shared Linux cluster |
| `standard` | Configurable | Repository default is SLURM/Singularity | Use only after inspecting its effective configuration |

Prefer an explicit named profile in generated commands. `setup.sh` can replace the `standard` block, so its behavior may differ between checkouts.

When Nextflow is installed, inspect rather than infer:

```bash
nextflow -version
nextflow config -profile <PROFILE>
```

Confirm Java and the chosen container runtime independently. A profile existing in `nextflow.config` does not prove that its runtime, scheduler, mounts, or hardware are available.

## Understand setup side effects

`setup.sh` is an installation/configuration action, not a harmless validation command. It can:

- validate an existing licensed ANNOVAR installation;
- download hg38 ANNOVAR databases;
- download and extract PhenoSV resources;
- download Phen2Gene knowledge-base assets;
- rewrite paths in downloaded resources;
- update ANNOVAR/PhenoSV host paths in `nextflow.config`; and
- replace `profiles.standard` with the selected backend.

Repeated setup can redownload resources. Run it only when the user asked to prepare or reconfigure PipeVar, and preserve intentional local configuration first. It does not install Nextflow, Java, Docker, Singularity, SLURM, ANNOVAR itself, a reference bundle, PhenoGPT2 models, or optional HmtVar data.

## Container mounts and cluster visibility

Configured host resources map to these container locations:

| Host parameter | Container destination |
| --- | --- |
| `--annovar_host_path` | `/annovar` |
| `--phenosv_host_path` | `/PhenoSV/train_data` |

The host paths must exist on the execution host. With SLURM, data, reference, model, output, work, and resource paths must be visible at consistent absolute locations on every eligible compute node.

When the resources already exist and only this run needs different mounts, the documented parameters can override the configured host paths without running setup:

```text
--annovar_host_path <ABSOLUTE_ANNOVAR_DIR>
--phenosv_host_path <ABSOLUTE_PHENOSV_DIR>
```

Validate that the selected profile uses those effective parameter values before launch. Treat a persistent configuration repair or resource download as a separately authorized setup task.

Do not expose full mount paths in public logs when they contain usernames, case identifiers, or restricted storage topology.

## PhenoGPT2 and GPU mode

`--GPU yes` changes supported GPU-aware processes, not only phenotype extraction. For short-read DeepVariant it selects the configured GPU image and runtime options; for PhenoGPT2 it supplies the required GPU container/scheduler settings. Confirm the selected profile sets the matching Docker or Singularity GPU backend and that the host actually provides a usable GPU. If neither DeepVariant nor PhenoGPT2 will run, omit `--GPU yes`; do not present it as general pipeline acceleration.

PhenoTagger is the default for clinical notes. Enable PhenoGPT2 only when explicitly requested and all required resources are configured:

```text
--phenotype_extractor phenogpt2
--GPU yes
--phenogpt2_model_host_path <ABSOLUTE_VERSIONED_NEW_MODEL_DIR>
```

With negation processing, also require:

```text
--phenogpt2_negation yes
--phenogpt2_negation_model_host_path <ABSOLUTE_NEGATION_MODEL_DIR>
--phenogpt2_embedding_model_host_path <ABSOLUTE_EMBEDDING_MODEL_DIR>
```

An optional pre-created cache uses `--phenogpt2_cache_host_path`. Model directories are mounted read-only; the cache is writable. Docker adds GPU runtime options and Singularity uses NVIDIA passthrough. SLURM GPU allocation comes from `--gpu_cluster_options`.

Current validation requires `--phenogpt2_batch_size 1`, `--phenogpt2_chunk_batch_size 1`, `--phenogpt2_max_forks 1`, and `--phenogpt2_wc 0`. Confirm future changes in [PARAMETERS.md](../../../../docs/PARAMETERS.md) and `main.nf` rather than changing them speculatively. HPO-only inputs do not require PhenoGPT2 models, cache, or GPU.

Use a fresh work directory whenever externally mounted model contents change. Nextflow cannot reliably detect mutations inside external mounts for cache validity.

## Resources and concurrency

PipeVar configures task-level CPU, memory, and time requests, with as many as three retries. Several processes increase CPU and memory on later attempts while keeping fixed time limits.

- Do not sum every process request into a whole-run minimum.
- Local capacity must cover the largest tasks that may execute concurrently.
- Use `--deepvariant_max_forks` to limit concurrent DeepVariant tasks when necessary.
- PhenoGPT2 is configured conservatively for concurrency; confirm its current fork setting before changing it.
- Use [BENCHMARKS.md](../../../../docs/BENCHMARKS.md) for planning context, not as a guaranteed runtime or capacity contract.

If a task exceeds local capacity, change execution capacity or deliberate concurrency/resource configuration. Do not blindly rely on retries, which can request more resources.

## Decide whether `-resume` is safe

Resume only when all of these remain unchanged:

- PipeVar code and included scripts/modules;
- input file paths and contents;
- pipeline parameter values;
- container identities/content;
- reference and database contents;
- externally mounted ANNOVAR, PhenoSV, PhenoGPT2, and other resource contents; and
- the original Nextflow work directory and cache metadata.

An infrastructure interruption, scheduler preemption, or transient runtime failure can be resumed after confirming the run contract is unchanged. A corrected input, changed reference/model/database, changed route/threshold, changed code, or removed `work/` is not the same cached run.

An existing nonempty output directory can be reused for a verified unchanged `-resume` when it is the original run's publish destination. Confirm that identity explicitly; otherwise treat nonempty output as a collision risk.

Resource-only scheduler adjustments can affect cache and execution behavior in ways that depend on the directive and Nextflow version. Follow the repository's strict reproducibility rule: treat changed configuration as a changed run unless the operator explicitly accepts and verifies cache reuse.

Before resuming:

1. diagnose the first failure with [troubleshooting.md](troubleshooting.md);
2. preserve `work/`;
3. reproduce the complete original command;
4. add `-resume` without silently changing other arguments; and
5. record why cache reuse is considered valid.
