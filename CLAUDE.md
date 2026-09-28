# Working in this repository

## Comments

One line inside the code body, if at all possible. Subroutine and module headers
may be longer but should still be as short as possible — what it does, and the
units and conventions a caller needs.

Never narrate history in a comment: no "previously this...", no "the sign was
flipped in 2016", no numbers from a debugging session. What was fixed is not what
the code does. That belongs in the commit message or in
`docs/assets/changes_fromV3Summa.txt`.

This applies to comments you add or touch. Leave pre-existing comments alone
unless asked.

## `docs/assets/changes_fromV3Summa.txt`

Every code change gets a numbered entry:

- **1–2 lines**, as short as possible. Match the entries numbered below 73 —
  those are the house style.
- Include a **short commit hash and a date**.
- Prefix the number with `-` when the change alters a backward-euler solution
  (i.e. is not backwards compatible).

```
 78) run_oneGRU.f90 cascade order built once per GRU, not every step, commit e9c746a3 Sep 25, 2026
-79) var_derive.f90 root density normalisation only trims an excess, commit 40e9929d Sep 24, 2026
```

Commit the code first so the hash exists, then add the entry, then commit the log.

## Git

- **Never** add a `Co-Authored-By:` line.
- Keep commit messages succinct: a clear subject, and a body only if something
  needs saying — a couple of lines, not paragraphs. No verification logs.
- **Never push.** Commit locally; the user pushes.

## Adding an `iLook` variable

`iLookFLUX`/`iLookPROG`/`iLookDIAG`/etc. in `build/source/dshare/var_lookup.f90`
are parameters built from a positional constructor of literal indices.

1. Append the member at the **end** of the type — never insert mid-type.
2. Extend the constructor with the next integer.
3. Register it in `dshare/popMetadat.f90`, `dshare/get_ixname.f90`, and
   `dshare/fluxMapping.f90` for fluxes (`state1=integerMissing,
   state2=integerMissing` for a pure diagnostic).

Skipping step 2 leaves the last member at its `integerMissing` default, which is
then used as an array subscript: SIGBUS during initialisation, no output.

## Building

An existing build directory needs these in the environment when it reconfigures:

```
FC=/opt/local/bin/gfortran LIBRARY_LINKS='-llapack' make -j8
```

## Tests

`utils/test/test_mflow/` holds five coupled Sagehen cases. Run them all after any
change to the coupling or the solver, and check the reported coupled water budget
against the previous run — `sagehen1` should report `-11971448.016178789` sent.
