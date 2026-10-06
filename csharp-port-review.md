# ALEXI Fortran → C# port — completeness review

**Date:** 2026-10-05
**Scope:** `csharp/` (Alexi.Core, Alexi.IO, Alexi.Cli, Alexi.Core.Tests, Alexi.Golden.Tests) vs `info/recipe/ALEXI_0.1/src/`
**Related docs:** [`csharp-conversion-plan.md`](csharp-conversion-plan.md) (plan of truth), [`handoff.md`](handoff.md) (execution state — **stale**, see §4.8)

---

## 1. Verdict

**The port is code-complete for `ALEXI_0.1`. It is not validated.**

The remaining work is almost entirely *proof*, not *translation*. No routine is missing; no numeric
decision has been checked against the original model.

---

## 2. Evidence (verified directly, not taken from `handoff.md`)

| Check | Result |
|---|---|
| `dotnet build Alexi.sln` | **0 errors, 0 warnings** (`TreatWarningsAsErrors=true`) |
| `dotnet test Alexi.sln` | **54 passed** — 49 `Alexi.Core.Tests` + 5 `Alexi.Golden.Tests` |
| Bare `dotnet build` / `dotnet test` | **ambiguous** — both `Alexi.sln` and `Alexi.slnx` exist |
| Model constants vs `ALEXI_parm.inc`, `USflux_parm.inc`, `flag.inc` | **exact match** (`mli=41`, `ml=mt=8000`, `kz=30`, `ilg=1456`, `jlg=625`, `kx=1440`, `ky=600`, `kt=8`, `bad=-9999`) |
| `GridDims.DirectRecord` vs `irec=(jj-1)*ilg+i`, `jj=jlg-j+1` | **exact match** — row flip preserved |
| `GetSoilHeat` / `GetSoilHeatOld` | **verbatim**, including the dead `g = 0.31*rnsoil` TESTING overwrite |
| `Program.cs` vs `USflux_main.f` | **exact match** — `call landcover(...)` → `call USflux(...)`, 4-arg contract |
| `real*4` fidelity | **clean** — 146 `MathF` calls, zero `double` declarations in `csharp/src` |
| CI workflow | **absent** — no `.github/` directory |
| Golden corpus | **absent** — no `tests/golden/`, no expected outputs, no sample inputs |

### 2.1 Call-graph fidelity is exact

Enumerated every `public`/`internal` method in `csharp/src/` that is **never called from `src/`** —
21 methods. Every one is dead or commented-out in the Fortran as well:

| C# method (unreachable) | Fortran status in `ALEXI_0.1` |
|---|---|
| `Cover.LoadTables` | `load_tables` never called (only `loadclasstable` is, from `landcover.f:25`) |
| `AlexiUSFluxClear.DailyFluxClearEF` | `daily_flux_clear_EF` never called; `clear_day_proc` calls `…_SDN` only |
| `AlexiInput.OpenOutput` | `open_output` never called |
| `BinaryGrid.ReadFlipped` / `ReadUnflipped` | `binread` / `binread_swap` never called |
| `AlexiCore.FindTaParallel` | `findta_parallel` never called |
| `AlexiCore.GetSoilHeatOld` | `getsoilheat_old` never called |
| `AlexiCore.GrowPbl2` | `growPBL2` never called |
| `AlexiInput.ExtractHourlyInputNLDAS` | commented out — `USflux_run.f:163` |
| `AlexiUSFluxUtl.CheckUSInput` | commented out — `USflux_run.f:306` |
| `AlexiInput.StoreOutputBin` | commented out — `USflux.f90:136`; C# gates it behind `--store-binary` |
| `AlexiOffline.PblRead` / `SfcRead` / `ArrayNav` / `WritePoint` | separate Fortran programs with their own `main` — no C# CLI entry (see §4.6) |

**Conclusion:** the C# mirrors the Fortran call graph *including its dead code*. These are faithful
dead routines, not missing wiring.

---

## 3. Remaining bugs

Ordered by severity. Items 1–4 are new findings not present in `handoff.md`.

### 3.1 Missing `USflux_main.f` stdout banner — **blocks golden stdout parity**

The `#### Begin Processing … ALEXI USA Version ####` block and the closing
`--->>  Processing is completed.` line are emitted nowhere in `csharp/src/`.
A golden test that diffs full stdout will fail on this before it fails on physics.

Also unverified: list-directed integer spacing. Fortran `write(6,*)'MDATE = ',MDATE` vs
`AlexiDriver.cs:27` `stdout.WriteLine("MDATE =  " + MDATE)`.

### 3.2 `Cover.LoadClassTable` — `tabdesc` is wrong twice (`Cover.cs:62`)

- Fortran **never reads** the file's `Class_description` column; it hardcodes `tabdesc(1..14)`
  MODIS biome names in `USflux_cover.f`. The C# takes the file's last column instead.
- The fallback `tok.Length > 14 ? string.Join(' ', tok.Skip(14)) : tok[13]` assigns **`fcmin`**
  as the description when a row has exactly 14 fields — off-by-one.
- **Impact today: zero.** `tabdesc` is declared in every Fortran `COMMON` but never appears in a
  `write`. Cosmetic — but it will bite any diagnostic output.

Column mapping for the other 14 fields was checked against `USflux_cover.f:57-62` and is **correct**,
including `itabclass(iclass)` ← column 12 (the MODIS biome class), which is why the C#'s initial
`itabclass[c] = iclass` is harmless dead code.

### 3.3 Culture-sensitive float formatting in `AlexiOffline.cs`

`W()`, `Join()` and the `us.input` block use bare `ToString()` / interpolated strings instead of the
`InvariantCulture` helpers used elsewhere. Currently masked only by `InvariantGlobalization=true` in
`csharp/Directory.Build.props` — a **build-level** mitigation, not a code-level one. `Alexi.IO` is a
library; any host without that flag emits comma decimal separators.

### 3.4 `GrowPbl` index clamp

Fortran: `it = int(dta/dtheta)`, guarded only for `it > 8000` → out-of-bounds / negative read for
small or negative ΔT. C# clamps `it < 1 → 1`. Documented in-code, but it is a real numeric
divergence on a reachable path and needs the author's decision.

### 3.5 `nclear` counter semantics

Fortran's `save` makes `nclear` sticky across `USflux` calls (and it is never initialized or
printed); C# resets it per call. Only matters for multi-day runs.

### 3.6 `pbl_read` / `sfc_read` ported but unreachable

`AlexiOffline.PblRead`, `SfcRead`, `ArrayNav`, `WritePoint` exist but the CLI exposes no subcommand
(only `interp-lai`). They also require a `GridState`, which eagerly allocates **684 MB**:

| Array | Size |
|---|---|
| `htht` (41×1456×625) | 149.2 MB |
| `ztan`, `zzan` (30×1456×625 each) | 109.2 MB each |
| `ctloc`, `cta`, `cea`, `cwind`, `csdn`, `cxlwdn`, `cpres`, `clst` (1440×600×8 each) | 27.6 MB each |
| **Total across `GridState` + `PixelState` + `DailyState`** | **683.6 MB** |

The CLI path never constructs `GridState`, so this cost is latent — but it must be made lazy or
pooled before these utilities are wired in.

### 3.7 `landcover.txt` data bug — faithfully reproduced

`info/recipe/ALEXI_0.1/store/landcover.txt` has 8 rows for classes 1–7 with class **5 duplicated**
(Cropland, then Wetland). Both Fortran and C# end up with class 5 = **Wetland**, not Cropland.
The port is correct; the input data is suspect. Confirm with the author.

### 3.8 `handoff.md` is stale

- Claims "nothing is committed / first commit pending" — three commits now exist (`42637f6`, `ac4756b`, `700d054`).
- Claims "unported routines throw `NotImplementedException`" — there are **zero** occurrences in `csharp/src`.

---

## 4. Test coverage gaps

All 54 existing tests are **hand-derived invariants and analytic checks**, not Fortran-parity
assertions. `Alexi.Golden.Tests` contains only `GetDate`/`GetDgmt` tests plus an empty
`Parity_Harness_Placeholder`.

Untested directly (~2,300 LOC of core physics and driver code):

- `AlexiCore.Alexi`, `Findhn`, `FindTa`, `FindTaParallel`, `GrowPbl`, `GrowPbl2`, `GetSoilHeat`
- `AlexiWater.AlexiWaterRun`
- `Cover.CoverProps`, `Cover.LoadTables`, `Cover.LoadClassTable`
- `AlexiUSFluxClear.ClearDayProc`, `DailyFluxClearSDN`, `DailyFluxClearEF`
- `AlexiInput.ExtractInput`, `StoreOutputBin`, `ScreenOutput`, `OpenOutput`
- `AlexiDriver.Landcover`, `AlexiDriver.USflux`
- `BinaryGrid` read/write round-trip including the row flip
- `InterpLAI`

---

## 5. Refactoring plan

### P0 — build the parity gate (this *is* the remaining task)

1. **Get a runnable Fortran oracle.** WSL2 is available on this host (`Ubuntu`, kernel 6.18.40.1);
   `bin/alexi_proc` is a Linux ELF (GNU/Linux 2.6.18, `with debug_info, not stripped`) and may run
   as-is. Otherwise `apt install gfortran scons` and rebuild with `-ffixed-line-length-132 -g`.
   Note: neither `gfortran` nor `scons` is currently installed in the WSL image.
2. **Create `tests/golden/`** — 3–5 cases (clear-sky, cloudy, water/ice, non-converging,
   bad-input pixel) with checked-in `met_` / `sat_` / `nparm_` / `profile_` / `veg_us_*.input`
   inputs and captured stdout + `.LOG` + `*_YYYYDDD.dat` grids + per-cell numeric fingerprint.
3. **Replace the placeholder** in `Alexi.Golden.Tests` with a real harness (Fortran as subprocess,
   C# in-process), tolerance `1e-3` absolute on fluxes.
4. **Add a CI workflow** — build + `dotnet test`. There is currently no `.github/` at all.
5. **Delete `Alexi.sln` or `Alexi.slnx`** so bare `dotnet build` / `dotnet test` is unambiguous.

### P1 — direct tests for the untested core

Cover §4 above, prioritising `AlexiCore.Alexi` → `Findhn` → `FindTa` → `GrowPbl` (the solver loop)
and `BinaryGrid` round-trip.

### P2 — de-duplicate the Fortran-formatting helpers

`F()`, `Pad()`, `Spaces()` are copy-pasted in `AlexiDriver.cs:344-353` and `AlexiInput.cs:439-451`
(plus `PadEnd`, `AlexiUtl.F92`, and `Pow`/`Cos` shims in `AlexiOffline.cs:225-227`). Extract one
internal `FortranFormat` class — this is also where the §3.3 culture fix belongs.

### P3 — only after parity is green

- `Parallel.For` over pixels
- `Span<float>` / `ArrayPool<float>` for the grids; lazy `GridState` allocation
- Structured logging; NuGet packaging — **resolve the BSD/MIT/CC0 license conflict first**

---

## 6. Still unported by design

- Python physics re-implementation (`ALEXI.py`, `ALEXI_utl.py`) — plan module 14, not started
- `pyalexi` iteration flags: `--start-year/--end-doy`, `--point`, `--region`
- `gzip`/`gunzip` shell-out in `pbl_read.f` / `sfc_read.f` — C# reads raw `.bin` only

---

## 7. Do not trust yet

- Green unit tests ≠ parity. Nothing has been compared against Fortran output.
- `FindAlbedoSoil` / `FindAlbedoVeg` convergence loops and `AlexiUtl.RunInit` are unvalidated.
- `DailyFluxCloudy` sets all daily fluxes to `BAD` — mirrors the source's gap-fill behaviour, not a
  computed cloudy flux; confirm it is intended.
- The `lookup` COMMON alias decision (`lookup_i`/`lookup_j` not reproduced) is unconfirmed; any
  golden case touching it may encode undefined behaviour.
- `interp_LAI.f` is **not in the SCons build** — decide whether it belongs in the port at all.
