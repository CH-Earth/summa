# Troubleshooting and preflight checks

Use this page to separate a configuration problem from an optimization result.
Do not remove constraints or suppress errors merely to make a job complete.

| Symptom | First checks |
| --- | --- |
| `{experiment}` is unresolved | Launch the manifest that defines `multi_case.experiment`, or use literal paths in a direct configuration |
| Basin inputs are missing | Compare `case_name` with the actual domain directory; check prefixes such as `CAN_`; inspect expanded paths |
| Observations are missing despite a zero target weight | Diagnostics still require valid observations and alignment; check explicit target `obs_file` and inherited path |
| An old example fails to parse | Replace `param_list`, separate bounds/transform tables, `sample_gaps`, and bare target selectors with the modern schema |
| Ratio field lengths differ | Check horizontal/domain correspondence, scalar versus layered native types, and whether explicit local `min` reduction is intended |
| Fixed-only constraint fails | Examine initialized physical fields; resampling unrelated variables cannot repair this relationship |
| Proposal retries are exhausted | Check incompatible bounds, minimum gaps, local minima, and a very small feasible volume |
| Many trials fail during initialization | Inspect the actual repeated error; do not assume every initialization failure is a bad parameter sample |
| Calibration output exists but results are missing | Check completion timestamps and `trial_failed`; allocated sample length alone does not prove execution |
| DDS diagnostic is missing but KGE is valid | An undefined zero-weight metric can legitimately be fill without failing the trial |
| NSGA-II rejects a diagnostic target | Zero-weight targets are currently DDS-only |
| A copied parameter vector gives different results | Distinguish coefficients from physical exports, and retain the same initial state, inputs, and observation processing |
| A multi-node direct run cannot read the shared restart | Only its filename is broadcast; the rank-zero restart must be reachable from every node |
| The job exceeds memory limits | Check full-series buffers, retained series, routing storage, and complete-instance replication across workers |
| More ranks do not improve throughput | Check trials per generation/budget, I/O contention, oversubscribed threads, and slow-trial tails |

## Before a production run

Confirm the actual executable and staged decisions file rather than inferring
them from a directory name. Check that forcing covers the automatic preceding
spin-up year and that observation periods match the intended evaluation dates.
Make sure parameter tables and trial inputs produce the intended initialized
fields, especially before spatial multiplication.

Verify the requested rank layout with a small job. Each group needs two ranks
at minimum, every trial currently holds the whole domain, and independent
repeats should use manifest run indices rather than repeated direct commands
with the same default seed.

Check that a completed case produces the expected target metadata, selected
physical export, and persistent copy. Retain the common state and original
configuration separately when they are needed for replay. A tiny trial budget
should not accidentally become the production budget, and a production budget
should not be used merely to diagnose a missing input path.

## During a run

Monitor completed evaluations and failures, not just process existence or CPU
allocation. Measure memory after spin-up and during full-period evaluations.
A transient early low-memory sample is not a full-run peak measurement.

Do not refresh or delete a node-local input/output tree while active workers are
using it. Inspect the error message associated with a failed trial before
attributing it to numerical instability. If all workers fail for the same input
reason, continuing thousands of proposals is not useful calibration.

## After a run

Read target dimensions and metadata before calculating maxima. Select DDS
parameters by the configured active objective and report diagnostics at that
sample. Keep failed/incomplete counts visible when comparing experiments.

For stochastic optimization, differences between independent repeats are part
of the result. A changed worker-completion order can change DDS proposals.
Use fixed-field replay for a numerical-equivalence question and matched basin,
period, budget, and repeat rules for a calibration-performance question.

See [Outputs](outputs.md) for interpretation and [Parallelization](parallelization.md)
for communicator and filesystem behavior.
