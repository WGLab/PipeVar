# Troubleshooting PipeVar runs

Diagnose the first real failure rather than the final cascade. Preserve the Nextflow `work/` directory and avoid changing parameters until the failing stage is understood.

## Capture the run contract

Record, with sensitive paths or identifiers redacted when shared externally:

- PipeVar Git revision and whether the checkout is dirty;
- complete launch command and working directory;
- selected profile and relevant effective configuration;
- Nextflow and Java versions;
- Docker or Singularity version and, for SLURM, scheduler context;
- reference build and reference/resource versions; and
- report/trace filenames and the first failing process name.

The configured Nextflow report and trace default to the launch directory; do not assume they are published under `--output_directory`.

Useful read-only commands include:

```bash
git status --short
git rev-parse HEAD
nextflow -version
nextflow config -profile <PROFILE>
```

Do not paste a whole clinical note, sample manifest, effective configuration, or unredacted absolute clinical-data paths into an issue.

## Find the first failure

1. Read `.nextflow.log` from the launch directory and locate the earliest process that failed, not merely a downstream channel cancellation.
2. Use the trace row or log entry to identify the task hash/work directory.
3. Inspect these files in that exact task directory:

   | File | Purpose |
   | --- | --- |
   | `.command.sh` | Exact generated task command |
   | `.command.run` | Nextflow task wrapper and staging behavior |
   | `.command.err` | Standard error from the task |
   | `.command.out` | Standard output from the task |
   | `.command.log` | Combined/runner log where present |
   | `.exitcode` | Task exit status where written |

4. Check staged input symlinks and their targets without modifying them.
5. Compare the failure to the route preflight and effective profile before deciding whether to rerun.

Avoid exposing PHI while quoting logs. Include the smallest excerpt that identifies the failure.

## Classify before acting

| Failure class | Typical evidence | Response |
| --- | --- | --- |
| Route validation | Top-level error before processes launch | Correct the unsupported or missing parameter combination; do not use resume as a substitute |
| Input/index/reference | Missing file, unreadable link, contig/build mismatch, decode/index error | Correct or regenerate the affected input only with authorization, then treat it as a changed run |
| Runtime/profile | Missing Docker/Singularity/SLURM command, local run submits to SLURM | Select or repair the intended explicit profile |
| Container | Pull failure, runtime permission, missing image | Verify runtime access and image identity; distinguish transient registry failure from invalid configuration |
| Bind/resource visibility | `/annovar`, `/PhenoSV/train_data`, model, or data path absent in task | Verify effective host path and node visibility; changing mounted content invalidates strict resume assumptions |
| Scheduler/resources | OOM, walltime, preemption, unavailable GPU | Separate transient preemption from inadequate requests; consider documented concurrency/resource controls |
| External models/databases | Missing/incompatible ANNOVAR, PhenoSV, Phen2Gene, or PhenoGPT2 files | Verify exact configured resource directory/version; do not run setup reflexively |
| Tool/process logic | Inputs are staged and runtime works, but caller/annotator exits deterministically | Inspect `.command.sh` and tool error, then reproduce only with non-sensitive minimal inputs where appropriate |
| Publishing/output | Task succeeds but expected artifact is absent or collision occurs | Check `publishDir`, selected route, output permissions, and conditional output contract |

## Common symptoms

### A local run tries to submit to SLURM

The default `standard` profile may use SLURM. Generate the command with `-profile local_docker` or `-profile local_singularity`, or intentionally reconfigure `standard` as a separate setup task.

### Containers cannot see ANNOVAR or PhenoSV

Inspect effective `annovar_host_path` and `phenosv_host_path`, then confirm the profile mounts them at `/annovar` and `/PhenoSV/train_data`. On a cluster, verify the same absolute host paths from compute nodes.

### A local task requests too many resources

Review the failing process and concurrent workload. `--deepvariant_max_forks` can limit DeepVariant concurrency. Later retry attempts may request more CPU/memory, so repeated automatic retries may worsen a capacity mismatch.

### PhenoGPT2 behaves inconsistently after model changes

Mounted model contents are external to normal input staging. Use a fresh work directory after changing model content and record the model version/path; do not rely on old cached tasks.

### Expected output is missing

First determine whether the selected route promises it. An HTML report, caller VCF, repeat table, or mitochondrial result can be conditional. Then confirm the relevant process completed and publishing targeted the expected `--output_directory`.

## Retry or resume decision

- Resume an unchanged run after confirmed scheduler preemption, transient runtime interruption, or another failure that did not require changing the run contract.
- Do not resume as though nothing changed after modifying input contents, reference/database/model contents, pipeline code, parameters, containers, or work/cache state.
- Do not delete `work/` while diagnosis or resume remains possible.
- Stop after the same deterministic failure repeats with unchanged evidence; continuing retries without a new hypothesis is not troubleshooting.

Use [runtime-and-resume.md](runtime-and-resume.md) for the complete cache-safety rule and [result-review.md](result-review.md) after a successful rerun.

## Reproducible handoff

A useful internal handoff contains:

1. failure classification and first failing process;
2. redacted command and relevant configuration values;
3. task hash/work path;
4. concise `.command.err`/`.command.out` excerpt;
5. input/reference/resource checks performed;
6. whether any run-contract element changed; and
7. recommended fresh-run or resume action with rationale.
