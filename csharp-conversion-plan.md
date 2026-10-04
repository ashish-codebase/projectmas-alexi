# Converting ALEXI / pyalexi to C#

**Target:** `projectmas_alexi` (ALEXI two-source ET model, Fortran + Python wrapper) → a .NET / C# solution.

**Measured scope** (`info/recipe/ALEXI_0.2/src`):

| Item | Count |
|---|---|
| Fortran source files (`.f`, `.f90`, `.inc`) | 26 |
| Fortran lines | 5,696 |
| Subroutines + 1 `program main` | 55 + 1 |
| `common /.../` declarations | 353 |
| `go to` statements | 26 |
| bare `save` statements (implicit static state) | 28 |
| `open(` / `close(` | 95 / 11 |
| `read(` / of those with `iostat` | 31 / **1** |
| `equivalence` sites | 4 (`USflux_run.f` ×2, `landcover.f` ×2) |
| `real*8` / `double precision` | **0** — model is entirely `real*4` |
| Python wrapper lines (`pyalexi.py`, `ALEXI.py`, `ALEXI_utl.py`) | ~700 |
| Automated tests | **0** (`info/test/run_test.py` is `import pyalexi`) |
| CI | none |

---

## 0. First decision: full port, or interop?

Three viable strategies. Pick before writing any C#.

| Strategy | What it is | Cost | When it wins |
|---|---|---|---|
| **A. Full port** | Rewrite all 5,700 Fortran lines in C# | High (~8–14 weeks solo, verified) | You need Windows-native, no gfortran, NuGet-distributable, parallel, testable |
| **B. Interop shell** | Keep `alexi_proc` (Fortran). Replace the Python wrapper with a C# CLI that shells out, or P/Invoke a Fortran shared library (`iso_c_binding`) | Low (~1–2 weeks) | You only need a C# front-end / .NET integration, not C# physics |
| **C. Transpile** | Run **Fortlang** (F#→C# converter) over the fixed-form source, then hand-fix residuals | Medium | You want C# source fast and can tolerate machine-generated idioms |

**Recommendation: A, executed as B→A.** Ship the interop shell first (it is cheap and gives you a working C# CLI on day 3), then migrate physics module-by-module behind the golden-output harness in Phase 0. This keeps a shippable artifact alive the whole time and never leaves you with a non-working model mid-port.

---

## 1. Phase 0 — Baseline and freeze (do this before any C#)

Without this you cannot tell a porting bug from a real change.

1. **Pick the source tree.** `ALEXI_0.1` and `ALEXI_0.2` are both in the repo. **DECIDED 2026-07: `ALEXI_0.1` is canonical** — `setup.py` (`version = "0.1"`) runs `scons` on `ALEXI_0.1/src`, `build.sh` runs `setup.py install`, so the shipped conda `bin/alexi_proc` is built from 0.1. Quarantine 0.2. Whitespace-insensitive diff shows the `.inc` files differ only in trailing whitespace (ported state/constants are valid for both trees), and `ALEXI_utl.f` is identical. Substantive `ALEXI.f` diffs (port the 0.1 version): 0.2 adds `rsmin` to `common/canopy` in `ALEXI.f`, adds `pres/ea` args to `findta` calls, comments out the `isnan(ta1)` stop (0.1 has it active), widens the iteration guard to `it.lt.0.or.it.gt.8000` (0.1: `it.gt.8000` only), and adds a VPD-based PT correction (`xlec = rndiv*pt*fg*(s/(s+0.66))`; 0.1 hardcodes `1.3`). `USflux_main.f`: 0.2 adds `character*3 arg3`.
2. **Build the reference binary** with gfortran (`-ffixed-line-length-132 -g`) on Linux, exactly as `SConstruct.py` does.
3. **Create a golden corpus:** 3–5 small input sets (`met_`, `sat_`, `nparm_`, `profile_`, `veg_us_*.input`) covering: clear-sky, cloudy, water/ice, non-converging, and bad-input pixels.
4. **Capture ground truth** for each case: stdout, the `.LOG` file, the text output on unit 103, and every `*_YYYYDDD.dat` direct-access binary.
5. **Record the numeric fingerprint** of each binary grid (per-cell `float`, checksum, max/min, count of `BAD = -9999`).
6. **Commit the corpus + expected outputs** to `tests/golden/`. These are the acceptance tests for every C# module.
7. **Fix the toolchain first:** remove committed artifacts (`*.o`, `.sconsign.dblite`, `__pycache__`, `.idea/`, `bin/alexi_proc`, `lib/`, `info/` build tree) from git; replace the `source: path: /PycharmProjects/projectMASalexi` absolute path in `meta.yaml`.

---

## 2. Phase 1 — Fix the defects that will not port

These are latent in Fortran (silent) and become loud in C# (exceptions). Fix in Fortran first, re-baseline, then port.

| # | Defect | Location | Why it must be fixed first |
|---|---|---|---|
| 1 | **Mismatched named COMMON `lookup`** — three incompatible layouts: `alv/aln/all/adv…` (`landcover.f`), `lookup_i/lookup_j` (`pbl_read.f`, `sfc_read.f`, `USflux.inc:27`, `USflux_run.f:48`), and `itabclass/tablai/…/tabdesc` (`USflux.inc:31`, `USflux_cover.f`, `USflux_run.f:435`) | `USflux.inc:27` vs `:31` | Undefined behaviour; the 7.28 MB `lookup_i/lookup_j` aliases the landcover tables. C# has no equivalent — you must decide the *intended* layout. |
| 2 | **`tabdesc` out-of-bounds write** — `tabdesc(nclass)` with `nclass=8`, written `tabdesc(1..14)` | `USflux_cover.f:57-70` | 300 bytes past the array. In C#: `IndexOutOfRangeException`. |
| 3 | **Landcover class hardcoded** — `iclass=5 !pvars(5)` | `USflux_run.f:97` | Every pixel runs as Cropland. Decide the intended behaviour before porting. |
| 4 | **`iflag` is two different entities** — `USflux.f90:14` declares `dimension iflag(ilg,jlg)` + `save`; `USflux_run.f:86,101,189,307` assign a scalar implicit integer local | `USflux.f90:86,107` | `iflag(ia,ja).eq.3` tests an **uninitialized** array. Water/ice and cloudy routing are unreliable. |
| 5 | **`ia`/`ja` uninitialized** — used in `call cover_props(ia,ja)` (`USflux_run.f:118`) and `writepoint`; only `ai`/`aj` are in `common/cover3` | `USflux_run.f:118` | C# requires definite assignment — will not compile. |
| 6 | **`screen_output` called with 0 args** against a 2-arg signature | `USflux.f90` → `USflux_run.f:407` | Interface mismatch; C# will not compile. |
| 7 | **`tabfcmin` all `-9999`** → `fcmin=-9999`, `if (fc.lt.fcmin)` never true; `perennial=.FALSE.` immediately overwrites the branch that sets it | `USflux_cover.f:141,159` | Dead code paths — port the *live* path only. |
| 8 | **`read(400,*)` with no sync check** — the `do while (iin.ne.ia…)` guard is commented out | `USflux_cover.f:187-190` | Silent misalignment of vegetation properties per pixel. |
| 9 | **No `iostat` on 30 of 31 reads, no error handling on 95 opens** | repo-wide | C# needs explicit exceptions; you must define the intended failure semantics. |
| 10 | **Hardcoded `/data/...` paths** (7 sites) and `./ALEXI/…` relative dirs | `USflux_run.f:768-770`, `sfc_read.f:85`, `landcover.f:311`, `USflux_dir.inc` | Must become injected configuration before porting. |
| 11 | **Unreachable code** after `go to 1000`; empty `if (writeme) then … endif` | `USflux.f90` | Delete so the port mirrors live behaviour. |
| 12 | **`interp_LAI.f` is in `Makefile` but not in `SConstruct.py`** | build files | Decide: port it, or drop it. |

Also: `ALEXI.py` and `ALEXI_utl.py` are currently **unimportable** (module-level call to undefined `profile`/`th2` at `ALEXI.py:196`; `Ts`/`Tc` undefined at `ALEXI.py:91`; `np`/`integrate` never imported in `ALEXI_utl.py`). `pyalexi.pyalexi:main()` references `args.start_doy`/`end_doy`/`start_year`/`end_year`/`point`/`region` — **none are defined in `arg_parse()`**, and the conda entry point calls `alexi()` with zero arguments. The Python layer is not a working reference; treat the Fortran binary as the only oracle.

---

## 3. Phase 2 — .NET scaffolding

```text
Alexi.sln
  src/
    Alexi.Core/            # physics, no I/O, no CLI
      Alexi.Core.csproj    # net8.0, NuGet package
    Alexi.IO/              # direct-access grids, list-directed readers
    Alexi.Cli/             # System.CommandLine front end
  tests/
    Alexi.Core.Tests/      # xUnit
    Alexi.Golden.Tests/    # parity vs Fortran corpus
  Directory.Build.props    # <TreatWarningsAsErrors>true</TreatWarningsAsErrors>
```

- **`net8.0`** (LTS) — needed for `System.IO.RandomAccess` and `System.Numerics.Tensors`.
- `Nullable enable`, `TreatWarningsAsErrors` — non-negotiable, because the Fortran's silent failure modes become compile-time signals here.
- NuGet: `System.CommandLine`, `xunit`, `FluentAssertions`, `MathNet.Numerics` (only if you port `ALEXI.py`'s numpy interpolation).
- Add a CI workflow (GitHub Actions): build + xUnit + golden parity. There is currently no CI at all.

---

## 4. Phase 3 — Translation rules (write this down, apply mechanically)

| Fortran in this repo | C# mapping | Trap |
|---|---|---|
| Implicit typing (`i–n`→integer, else real), no `implicit none` | Explicit `int`/`float` on every symbol | This is where most hidden bugs surface. Do it file-by-file. |
| `real*4` (the *only* real kind here) | `float` (IEEE binary32) | **Do not promote to `double`.** It changes results and breaks parity. Keep `float` end-to-end. |
| `integer*8 filg, fjlg` (`flag.inc`) | `long`, or `int` constants (1456, 625 fit in `int`) | Only place with 8-byte ints. |
| `common /name/ a,b,c` (353 decls) | One `AlexiState` class, grouped by block: `State.Cover`, `State.Data1`, `State.Data2`, `State.Pbl`, `State.Flags`… | Pass `AlexiState` explicitly to every routine. **No `static` mutable fields** — that is how you keep the model testable and parallel-safe. |
| `include 'X.inc'` (textual) | `Constants` / `GridDims` static classes + the `AlexiState` partials | `.inc` files carry `parameter`, `data`, and `common` — split them, don't transliterate. |
| `save` on subroutine locals (28 sites) | Instance fields on a per-run `AlexiState` | `save` with no entity list saves *everything* → today the model is single-threaded by construction. |
| 1-based arrays, e.g. `htht(41,1456,625)` | 0-based `float[]` (flattened) or `Span<float>` | Every index needs a `-1`. Flatten to 1-D: `idx = (j*ILG) + i`. Prefer flat arrays — 440 MB of multi-dim `common` arrays must be heap-allocated, not on the stack. |
| `go to` / labels (26) | `while`/`break`/labeled loops; C# `goto` is legal for the 3–4 genuinely irreducible cases | `USflux.f90`'s `go to 900` → re-test of `converged` is a control-flow bug, not a construct. Restructure. |
| `equivalence (s,ii)` (4 sites, byte-swap) | `BinaryPrimitives.ReverseEndianness` / `MemoryMarshal` | Delete `byteswapr4` entirely; use `BinaryPrimitives.WriteSingleLittleEndian`. |
| `open(unit=n, form='unformatted', recl=4, access='direct')` | `System.IO.RandomAccess` — `WriteSingleLittleEndian(handle, offset, value)` | Offset = `(record-1) * 4`. Record formula is `irec = (jlg-j)*ilg + i` (`USflux_run.f:920`) — note the **row flip** `jj = jlg-j+1`. Preserve it exactly. |
| `read(u,*)` list-directed | Hand-rolled `ListDirectedReader` (whitespace-delimited, comma = field boundary, EOF → exception) | Fortran list-directed has real semantics; a naive `Split()` diverges on edge cases. |
| `write(u,1003) 'label', value` with `format(a27,f12.3)` | `string.PadLeft/PadRight` + `value.ToString("F3", CultureInfo.InvariantCulture)` | **`CultureInfo.InvariantCulture` everywhere** — decimal comma locales corrupt output. |
| `stop 'msg'` (4) | Custom exceptions: `AlexiInputException`, `AlexiConvergenceException` | Never `Environment.Exit` from physics code. |
| `data` statements, `parameter` | `static readonly` fields / `const` | `USflux_dir.inc` `data HOMEDIR/"./ALEXI/"` → `IPathProvider`, not a constant. |
| Fortran character semantics (blank-padded, `index(s,' ')-1`) | `string` + explicit trimming | `l1=index(TABDIR,' ')-1` yields `-1` when the buffer is full → invalid substring. Replace with `TrimEnd()`. |
| `mod`, `aint`, `abs`, `alog`, `exp`, `cos` | `%`, `(float)Math.Truncate(...)`, `MathF.*` | Use `MathF` (float) to preserve precision, not `Math` (double). |
| `if (x.eq.-9999.)` on floats | `const float Bad = -9999f` + exact compare | Fortran does exact float compare here; C# `==` on `float` matches. Do **not** "fix" it with tolerance — that changes results. |

---

## 5. Phase 4 — Port order (bottom-up, one module per PR, parity-gated)

Each step: port → compile → run golden case → diff vs Fortran → merge.

1. **`Constants` + `GridDims` + `AlexiState`** ← `ALEXI_parm.inc`, `USflux_parm.inc`, `ALEXI.inc`/`ALEXI_f90.inc`, `USflux.inc`, `USflux_grids.inc`, `flag.inc`, `date.inc`, `USflux_dir.inc`. Resolve the `lookup` COMMON conflict here (Phase 1 #1).
2. **`AlexiUtl`** ← `ALEXI_utl.f`: `set_constants`, `checkvalue`, `ALEXI_errorcode`, `getprofile`, `getpbltable`, `runinit`, `updatefc`, `checkinput`.
3. **`Radiation`** ← `ALEXI_rad.f`: `getsunzen`, `getradprops`, `getnetrad`, `getnetradnight`, `getxlwdn`, `find_albedo_soil`, `find_albedo_veg`.
4. **`Atmos`** ← `ALEXI_atmos.f`: `getresistance`, `psimhn`, `psimhn_ky`, `canopyarch`.
5. **`Cover`** ← `USflux_cover.f`: `load_tables`, `cover_props`, `getfveg` (+ `store/landcover.txt` reader).
6. **`TwoSource`** ← `ALEXI.f`: `findhn`, `growPBL`, `growPBL2`, `findta`, `findta_parallel`, `getsoilheat`, `checksoln`. **Hardest module** — iterative convergence with `stopiter`/`converged` flags.
7. **`Water`** ← `ALEXI_water.f`: `alexi_water`, `findflux_water`, `findta_water`, `getsoilheat_water`, `getnetrad_water`.
8. **`DailyFlux`** ← `USflux_clear.f`, `USflux_cloud.f`, `USflux_rad.f`: `clear_day_proc`, `cloudy_day_proc`, `daily_flux_clear_SDN`, `daily_flux_clear_EF`, `daily_flux_cloudy`, `getradcomps`, `getnetrad_simple`.
9. **`Alexi.IO`** ← `USflux_run.f` I/O half: `binopen`, `binopen2`, `binwrite`, `binread`, `binread_swap`, `open_output`, `store_output_bin`, `screen_output`.
10. **`Input`** ← `USflux_run.f` input half: `extract_input`, `getdgmt`, `extract_hourly_input_nldas`.
11. **`Driver`** ← `USflux.f90` + `USflux_main.f`: argument parsing, the pixel loop, counters (`ntot/nbad/nconv/nfail/ncloud/nwater`), logging.
12. **Offline utilities** ← `pbl_read.f`, `sfc_read.f`, `landcover.f`, `interp_LAI.f` (decide scope; `interp_LAI.f` is not in the SCons build).
13. **CLI** ← `pyalexi/pyalexi.py` → `Alexi.Cli` with `System.CommandLine`. Fix the argparse defects while porting: define `--start_doy`/`--end_doy`/`--start_year`/`--end_year`/`--point`/`--region` (referenced but never declared); replace `nargs='*', type=bool` (which makes `bool("False") == True`) with proper `bool` flags; check the child process exit code.
14. **Python physics reimplementation** ← `ALEXI.py` + `ALEXI_utl.py` → C#. **Asset, not liability:** these are vectorized numpy re-implementations of the same two-source physics. Port them as an independent cross-check against the Fortran-derived `TwoSource`, and use the divergence to find porting errors. Fix their bugs first (`ALEXI.py:91` `Ts`/`Tc`; `ALEXI.py:196` module-level call; missing `numpy`/`scipy` imports in `ALEXI_utl.py`; `scipy` is undeclared in `meta.yaml`).

---

## 6. Phase 5 — Configuration and data layout

Replace every hardcoded path with injected config. Today: `./ALEXI_LOG/`, `./INPUTS/ALEXI_INPUT/`, `./ALEXI/store/landcover.txt`, `/data/data123/chain/4KM/GBIM/fluxes/`, `/data/data123/chain/4KM/GBIM/inputs/`, `/data/chain/CONUS/ALB/`, `/data/data123/chain/4KM/GBIM/inputs/MLAI`, `/data/data123/chain/4KM/GBIM/inputs/CLASS.dat`.

```csharp
public sealed record AlexiPaths(
    string InputDir, string OutputDir, string TableDir,
    string LogDir, string AlbedoDir, string LaiDir);

public sealed record AlexiRunOptions(
    int Mdate, string TimeTag, string Part, int NPoints,
    AlexiPaths Paths, int GridI = 1456, int GridJ = 625);
```

- Resolve paths relative to a config file / env var, never the CWD.
- Ship `store/landcover.txt` as an embedded resource **and** allow an override path.
- Note the repo's `store/landcover.txt` has 8 rows for classes 1–7 with class **5 duplicated** (Cropland, then Wetland) — the second row silently overwrites class 5. Encode the intended class table explicitly.

---

## 7. Phase 6 — Validation gate (the step that decides success)

1. **Golden parity tests** — for each corpus case, assert C# output equals Fortran output within a stated tolerance. Start at `1e-3` absolute on fluxes, tighten toward `1e-5` as you localize divergence. Expect bit-exact agreement on the pure-arithmetic modules; convergence-sensitive modules (`findhn`, `findta`, `growPBL`) will diverge first.
2. **Binary grid parity** — compare every `*_YYYYDDD.dat` cell-by-cell, including the `BAD = -9999` pattern and the `jj = jlg-j+1` row flip.
3. **Per-module parity** — expose internal state (`ta1`, `ts1`, `tc1`, `ra2`, `rs2`, `rx2`, `z2`, `hn`, `psi`) as a debug dump so you can bisect divergence to a single routine.
4. **Property tests** on the pure functions: `getfveg`, `psimhn`, `getsunzen`, `getxlwdn` — monotonicity, domain, and range invariants.
5. **Perf baseline** — time the Fortran run per pixel; require the C# version to be within 1.5× before parallelizing, then measure the speedup from `Parallel.For` over the pixel loop.
6. **Only after parity is green** may you refactor for design (dependency injection, `Span<T>`, `struct` state).

---

## 8. Phase 7 — Post-parity wins (only after green)

- **Parallel pixel loop:** `Parallel.For` over `(i,j)`. Requires Phase 1 #4/#5 fixed and **zero `static` mutable state** — the 28 bare `save` statements and 353 `common` blocks are what make the current code single-threaded.
- **`Span<float>` / flat arrays** for the ~440 MB of grid arrays (`htht(41,1456,625)` alone is 149 MB) — allocate once, rent from `ArrayPool<float>`, avoid LOH churn.
- **Structured logging** (`Microsoft.Extensions.Logging`) replacing unit-99 writes; close the log deterministically (the Fortran `close` calls are commented out at `USflux.f90:148-156`).
- **Real CLI contract:** `--start-year/--end-year/--start-doy/--end-doy` iteration (the dead `num_parts`/`parts`/`doys`/`years` computation in `pyalexi.py:main()` is the unfinished intent), `--version`, exit codes, `--dry-run`.
- **NuGet packaging** replaces conda + SCons + `setup.py`. Resolve the license conflict first: `setup.py` says BSD 3-Clause, classifiers say MIT, `meta.yaml` says `CCO` (typo of CC0) with family BSD. Pick one.

---

## 9. Tooling

- **Fortlang** (`f# → C#` transpiler) — the only serious automated option for fixed-form F77. Expect it to handle arithmetic, FORMAT, and I/O scaffolding; expect it to choke on the mismatched COMMON, `equivalence`, and `include`-based state. Use it for a head start on Phases 2–4, then hand-verify.
- **f2c → C → C++/CLI** — legacy, produces unreadable output. Not recommended.
- **Manual port with the golden harness** — recommended for the physics modules (`ALEXI.f`, `ALEXI_atmos.f`, `ALEXI_rad.f`).

---

## 10. Effort and risk

| Phase | Effort (solo, verified) |
|---|---|
| 0 — baseline + golden corpus | 1–2 weeks |
| 1 — Fortran defect fixes + re-baseline | 1–2 weeks |
| 2–3 — scaffolding + translation rules | 2–3 days |
| 4 — port modules 1–12 with parity | 6–10 weeks |
| 5–6 — config, CLI, packaging | 1 week |
| 7 — parallelism + cleanup | 1–2 weeks |
| **Total (Strategy A)** | **~10–16 weeks** |
| **Strategy B (interop shell only)** | **~1–2 weeks** |

### Top risks, ranked

1. **The `lookup` COMMON conflict** — you cannot transliterate it; you must decide intended semantics. Highest scientific risk.
2. **`iclass=5` hardcode** — the reference outputs may encode a debug hack as "truth". Confirm with the model author before treating it as the baseline.
3. **Convergence-path divergence** in `findhn`/`findta`/`growPBL` — single-precision iteration order is fragile; parity may be unreachable bit-exactly. Budget time for tolerance-setting.
4. **No tests, no CI, no reference inputs in the repo** — the golden corpus must be manufactured from scratch.
5. **Two source trees** (`0.1` built by `setup.py`, `0.2` by `SConstruct.py`) — confirm which is authoritative before baselining.

---

## 11. Minimum viable first sprint (concrete)

1. Declare `ALEXI_0.2` canonical; strip committed artifacts from git.
2. Build `alexi_proc` with gfortran; run it on one small input part; capture stdout + `.LOG` + all `.dat` grids.
3. Create `Alexi.sln` with `Alexi.Core`, `Alexi.Cli`, `Alexi.Golden.Tests`.
4. Port `Constants` + `GridDims` + `AlexiState` + `set_constants` + `checkvalue`.
5. Write the golden parity test harness (Fortran as subprocess, C# in-process).
6. Ship `Alexi.Cli` in interop mode (shells out to `alexi_proc`) — a working C# CLI on day 3.
7. Then start Phase 4 at module 2 (`AlexiUtl`).
