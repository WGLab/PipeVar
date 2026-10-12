# Nextflow compatibility

PipeVar requires Java 17 or newer and Nextflow 24.10.6 or newer. It supports
the legacy parser used by Nextflow 24.10/25.x and the native v2 parser used by
Nextflow 26.04. No parser environment variable is required for supported
default configurations.

## Verified matrix

The following matrix was verified on 2026-10-07 with Eclipse Temurin
17.0.20.1+1. Legacy lanes passed 22 harness checks and strict-parser lanes
passed the same checks plus recursive `nextflow lint` (23 checks).

| Nextflow | Build | Parser selection | Profiles | Lint | Contracts/smoke | Result |
| --- | ---: | --- | --- | --- | --- | --- |
| 24.10.6 | 5937 | Default (v1) | All four pass | Not available | Pass | Supported |
| 25.04.8 | 5956 | Default (v1) | All four pass | Not available | Pass | Supported |
| 25.10.7 | 12755 | Explicit v1 | All four pass | Not run | Pass | Supported |
| 25.10.7 | 12755 | Explicit v2 | All four pass | Pass | Pass | Supported |
| 26.04.6 | 12646 | Default/native v2 | All four pass | Pass | Pass | Supported |
| 26.04.6 | 12646 | Explicit v1 fallback | All four pass | Not run | Pass | Temporary fallback only |

The four profiles are `standard`, `slurm_singularity`, `local_singularity`,
and `local_docker`. The contract result covers structure, vendored nf-core
integrity, strict-syntax rules, compatibility-sensitive behavior, help
semantics, channel joins, and valid and invalid manifest boundaries. The smoke
process imports `INPUT_PREPARATION`, selects the single-sample VCF/SNP route,
executes one host-local task, and asserts `SMOKE|vcf|false|snp` without pulling
a container.

Verified one-JAR SHA-256 values:

| Nextflow | SHA-256 |
| --- | --- |
| 24.10.6 | `82c2a0690f830e6cee1fcb6fdd98d35c268f1d71a765444e796ab2643af3fcac` |
| 25.04.8 | `c73e44801457136872200c80e70152950a4663c11ba15ff58e6584719c52c099` |
| 25.10.7 | `929cb86438ddae07da2ae80f44c093be1ab6c216770a880ecafa22d534503610` |
| 26.04.6 | `2ca0251ae2d749317d9fbe5fe191a1616b5f44b608224268924c71b32f5ed9e2` |

## Reproduce the checks

Provide an official Nextflow launcher and a Java 17+ runtime, then run:

```bash
NEXTFLOW_BIN=/absolute/path/to/nextflow \
JAVA_HOME=/absolute/path/to/jdk-17 \
bash scripts/check_nextflow_compatibility.sh
```

The harness uses separate `NXF_HOME`, launch, work, and log directories for
each lane. By default it writes a retained result tree below `/tmp` and records
the PipeVar revision and source status, commands, exit codes, Java metadata,
Nextflow build, parser selection, runtime hashes, logs, and results. Override
the result location or lane matrix when needed:

```bash
COMPAT_ROOT=/absolute/writable/path \
NEXTFLOW_COMPATIBILITY_LANES='25.10.7@v2 26.04.6@default' \
NEXTFLOW_BIN=/absolute/path/to/nextflow \
JAVA_HOME=/absolute/path/to/jdk-17 \
bash scripts/check_nextflow_compatibility.sh
```

The launcher downloads missing official runtime artifacts unless its
environment is offline. For an offline run, populate each isolated
`NXF_HOME/framework/<version>` directory first.

## Interpreting failures

| Failure | Classification | Action |
| --- | --- | --- |
| Java reports class-file version 61 but supports only 55 | Host prerequisite | Use Java 17 or newer. |
| `.nextflow/history.lock` or bootstrap temporary file is unwritable | Host prerequisite | Use writable launch and `NXF_HOME` directories. |
| `sbatch` is missing | Profile readiness | Use a local profile or run on a configured SLURM host. |
| Docker socket access is denied | Container runtime | Correct Docker access or select a working Singularity profile. |
| ANNOVAR or PhenoSV mount directory is missing | Resource readiness | Run setup or correct the configured host path. |
| A local task requests more memory than the host provides | Host capacity | Use suitable compute or a test-only fixture override; do not lower production requests to satisfy compatibility tests. |
| An invalid manifest, phenotype combination, VCF mode, or mitochondrial route is rejected | Expected validation | Confirm the documented diagnostic; this is a passing negative test. |
| Nextflow 26.04 reports v2 parser syntax errors | Workflow compatibility | Treat it as a code defect and preserve the harness log. As a temporary diagnostic fallback, compare with `NXF_SYNTAX_PARSER=v1`; do not make v1 the normal v26 configuration. |
| A DSL2/module/config error remains with Java 17 | Workflow compatibility | Treat it as a code defect and preserve the harness log for investigation. |

This matrix does not run complete clinical routes, pull analysis containers,
compare biological results, or validate site-specific Docker, Singularity, or
SLURM installations.
