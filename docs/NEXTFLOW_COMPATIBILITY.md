# Nextflow compatibility

PipeVar requires Java 17 or newer and Nextflow 24.10.6 or newer. The workflow
currently uses the v1 DSL2 syntax parser. Set `NXF_SYNTAX_PARSER=v1` explicitly,
especially with Nextflow 26.04, where the v2 parser is the default.

## Verified matrix

The following matrix was verified on 2026-09-30 with Eclipse Temurin
17.0.16+8. Each release passed all 13 harness checks.

| Nextflow | Build | Profile configuration | Graph preview | Contract suites | Local smoke | Result |
| --- | ---: | --- | --- | --- | --- | --- |
| 24.10.6 | 5937 | All four pass | Pass | Pass | Pass | Supported |
| 25.04.8 | 5956 | All four pass | Pass | Pass | Pass | Supported |
| 25.10.7 | 12755 | All four pass | Pass | Pass | Pass | Supported |
| 26.04.6 | 12646 | All four pass | Pass | Pass | Pass | Supported with `NXF_SYNTAX_PARSER=v1` |

The four profiles are `standard`, `slurm_singularity`, `local_singularity`,
and `local_docker`. The contract result covers structure, vendored nf-core
integrity, compatibility-sensitive behavior, channel joins, and valid and
invalid manifest boundaries. The smoke process imports `INPUT_PREPARATION`,
selects the single-sample VCF/SNP route, executes one host-local task, and
asserts `SMOKE|vcf|single|snp` without pulling a container.

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
each release. By default it writes a retained result tree below `/tmp` and
records commands, exit codes, Java metadata, the selected syntax parser, and
runtime hashes. Override the result location or matrix when needed:

```bash
COMPAT_ROOT=/absolute/writable/path \
NEXTFLOW_VERSIONS='24.10.6 25.04.8 25.10.7 26.04.6' \
NEXTFLOW_SYNTAX_PARSER=v1 \
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
| Nextflow 26.04 reports v2 parser syntax errors | Parser mismatch | Export `NXF_SYNTAX_PARSER=v1`; native v2 migration is not part of the verified interface. |
| A DSL2/module/config error remains with Java 17 and parser v1 | Workflow compatibility | Treat it as a code defect and preserve the harness log for investigation. |

This matrix does not run complete clinical routes, pull analysis containers,
compare biological results, or validate site-specific Docker, Singularity, or
SLURM installations.
