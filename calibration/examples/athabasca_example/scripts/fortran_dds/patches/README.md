# Patches to SUMMA develop's Fortran DDS

Both patches change only `build/source/calibration/summa_parameter_sampling.f90` and are made against the
develop version checked out on Oct 23, 2026 and used for the original run (job 21621404). Apply one of them
in a checkout of that version:

    cd summa-develop
    patch -p1 < .../fortran_dds_seed.patch            # or fortran_dds_seed_defstart.patch

| Patch | Adds | Used for |
|---|---|---|
| `fortran_dds_seed.patch` | seed from `SUMMA_DDS_SEED` (default 42), printed as `DDS random seed = n` | seeds 43-45 and 43b-45b |
| `fortran_dds_seed_defstart.patch` | the seed option, plus `SUMMA_DDS_START=default`: the first sample is the SUMMA default set and DDS steps start from it; the 15 starting values are printed | default-start runs d42-d44 |

With neither variable set, both patched programs behave exactly like the unpatched code (seed 42, random start).
The default-start build was installed as a separate executable (`summa_mizuroute_opt_defstart.exe`);
see `slurm/run_fortran_dds_defstart.slurm`.

**Current develop:** the file was rewritten upstream on Nov 9, 2026 (merge of PR #674), so these patches no
longer apply to it. The new version already sets the seed through `run_index` (`seed = 42 + run_index - 1`),
which replaces `fortran_dds_seed.patch`; it has no default-start option yet.
