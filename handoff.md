# Handoff — ALEXI (Fortran) → C# port

**Date:** 2026-10-03 · **Branch state:** `csharp/` and `csharp-conversion-plan.md` are **untracked** in git (last commit `42637f6` is docs-only).
**Plan of truth:** [`csharp-conversion-plan.md`](csharp-conversion-plan.md) — this file records *execution state against that plan*, not a new plan.

---

## 1. Status snapshot (verified today)

| Check | Result |
|---|---|
| `dotnet build Alexi.slnx` / `Alexi.sln` (full solution) | ✅ **0 errors** (legacy `.sln` added, opens in VS Community) |
| `dotnet test` → `Alexi.Core.Tests` | ✅ **49 passed, 0 failed** |
| `dotnet test` → `Alexi.Golden.Tests` | ✅ **5 passed, 0 failed** (GetDate/GetDgmt; parity cases still empty) |
| CLI smoke run (synthetic 1-pixel inputs) | ✅ full pipeline ran: `veg_` → `veg_us` → cloudy path → `output_*` line + `.LOG` counters (no crash) |
| Golden corpus (`tests/golden/`) | ❌ does not exist |
| Fortran oracle | `bin/alexi_proc` is a **Linux ELF** — cannot be executed on this Windows host |

**Net:** the whole pipeline (`landcover` → `USflux` → `extract_input` → solver → daily flux → screen/binary output) is ported and builds; unit tests are green but there is **still no parity gate**.

LOC so far: ~5,056 lines across `src/` + `tests/` (24 source files, 8 test files).

---

## 2. What has been done

### Phase 0 — baseline (partially done, decisions locked)
- ✅ Canonical source tree **decided: `ALEXI_0.1`** (`setup.py` → `scons` → shipped `bin/alexi_proc`). `ALEXI_0.2` quarantined; substantive diffs catalogued in the plan.
- ❌ Golden corpus **not built** — no reference inputs, no captured stdout/`.LOG`/`*.dat`, no numeric fingerprint. This is the single largest remaining gap.

### Phase 1 — Fortran defects (catalogued, not fixed in Fortran)
- 12 defects enumerated in the plan; **handled at the C# boundary** rather than by patching Fortran:
  - `hourly_flux(ibad)` arity mismatch → documented in `AlexiUSFluxUtl.cs` header, extra arg ignored.
  - `lookup` COMMON conflict → **resolved by separation** (see §3 below).
  - `iclass` hardcode, `iflag` uninitialized array, `ia/ja` unassigned → surfaced as explicit state fields; not yet re-decided with the author.

### Phase 2 — .NET scaffolding (done)
- `csharp/Alexi.slnx` with `Alexi.Core`, `Alexi.IO`, `Alexi.Cli`, `Alexi.Core.Tests`, `Alexi.Golden.Tests`.
- `Directory.Build.props`: `net8.0`, `Nullable=enable`, `ImplicitUsings`, `TreatWarningsAsErrors=true`, `InvariantGlobalization=true`.

### Phase 3–4 — translation rules + module port order (plan §5 order)

| # | Module | Status |
|---|---|---|
| 1 | `GridDims` (`ALEXI_parm.inc`, `USflux_parm.inc`, `flag.inc`) + `PixelState`/`DailyState`/`GridState`/`LandcoverTables` + `AlexiPaths`/`AlexiRunOptions` | ✅ done — column-major flattening, `Bad = -9999f` exact, `DirectRecord` keeps the `jj = jlg-j+1` row flip |
| 2 | `AlexiUtl` ← `ALEXI_utl.f` (`SetConstants`, `GetProfile`, `GetPblTable`, `RunInit`, `UpdateFc`, `CheckInput`, `CheckValue` ×2, `AlexiErrorcode`) | ✅ done |
| 3 | `AlexiRad` ← `ALEXI_rad.f` (`GetSunzen`, `GetRadProps`, `GetNetRad`, `GetNetRadNight`, `GetXlwdn`, `FindAlbedoSoil`, `FindAlbedoVeg`) | ✅ done |
| 4 | `AlexiAtmos` ← `ALEXI_atmos.f` (`GetResistance`, `Psimhn`, `PsimhnKy`, `CanopyArch`) | ✅ done |
| 5 | `Cover` ← `USflux_cover.f` (`load_tables`, `cover_props`, `getfveg` 100-step grid search, `loadclasstable`) + `landcover.f` stage (`weighted_avg`, `canopyarch_lai`) | ✅ done |
| 6 | `TwoSource` ← `ALEXI.f` (`Alexi`, `Findhn`, `FindTa`, `FindTaParallel`, `GrowPbl`, `GrowPbl2`, `GetSoilHeat`, `CheckSoln`) | ✅ done |
| 7 | `Water` ← `ALEXI_water.f` (`AlexiWaterRun`, `FindFluxWater`, `FindTaWater`, `GetSoilHeatWater`, `GetNetRadWater`) | ✅ done (dead path in shipped 0.1: `iswater_inland` always false) |
| 8 | `DailyFlux` ← `USflux_clear.f` / `USflux_cloud.f` / `USflux_rad.f` | ✅ done (`ClearDayProc`, `DailyFluxClearSDN`, `DailyFluxClearEF` added) |
| — | `AlexiUSFluxUtl` ← `USflux_utl.f` (`HourlyFlux`, `Faopm`, `Integrate`, `GetTdepart2`, `CheckUSInput`, `SetPBLHeights`) | ✅ done |
| 9 | `Alexi.IO` ← `USflux_run.f` I/O half (`BinaryGrid` reader/writer, `OpenOutput`, `StoreOutputBin`, `ScreenOutput`) | ✅ done |
| 10 | `Input` ← `USflux_run.f` input half (`ExtractInput`, `GetDgmt`, `ExtractHourlyInputNLDAS`, `ListDirectedReader`) | ✅ done |
| 11 | `Driver` ← `USflux.f90` + `USflux_main.f` + `landcover` stage + `getdate` | ✅ done (`AlexiDriver.cs`, `RunCounters`) |
| 12 | Offline utils `pbl_read.f`, `sfc_read.f` (`AlexiOffline.cs`: `PblRead`, `SfcRead`, `ArrayNav`, `WritePoint`) | ✅ done (gunzip/gzip shell-out NOT reproduced — reads raw `.bin`) |
| 12b | `InterpLAI` ← `interp_LAI.f` | ✅ done (compile errors fixed) |
| 13 | CLI ← `pyalexi.py` | ⚠️ `USflux_main.f` 4-arg contract + dir overrides + `--store-binary` + `interp-lai` done; pyalexi iteration flags (`--start-year/--end-doy`, `--point`, `--region`) not implemented |
| 14 | Python physics re-impl (`ALEXI.py`, `ALEXI_utl.py`) | ❌ not started |

### Tests written (54, all green)
`GridDimsTests` · `AlexiUtlTests` · `AlexiRadTests` · `AlexiAtmosTests` · `AlexiUSFluxRadTests` · `AlexiUSFluxUtlTests` · `NewCoreTests` (GetFveg, CheckSoln, GetNetRadWater) · `DriverTests` in Golden suite (GetDate, GetDgmt).
These are **hand-derived invariants and analytic checks** — they are **not** Fortran-parity assertions.

---

## 3. Decisions locked in (do not re-litigate without new evidence)

- **Canonical tree = `ALEXI_0.1`** (overturns the plan's §11 "declare 0.2 canonical" line; `setup.py`/`build.sh` evidence is what decided it).
- **`real*4` end-to-end** — `float` + `MathF`, never `double`. Promotion would break parity.
- **No `static` mutable state** — every routine takes `PixelState`/`DailyState`/`GridState`/`LandcoverTables` explicitly. Keeps the model parallel-safe and testable; this is the deliberate answer to the 353 `common` blocks and 28 bare `save`s.
- **1-based Fortran indices preserved at call sites**, flattened via `GridDims.Index2/Index3*`; `DirectRecord` keeps the row flip verbatim.
- **`lookup` COMMON conflict resolved by separation**: `LandcoverTables` models only the `USflux.inc:31` layout (`itabclass/tablai/…/tabdesc`); the `lookup_i/lookup_j` aliasing (`USflux.inc:27`) is **not** reproduced. Documented in-file as unverified — needs author confirmation before any reference output that depends on the alias is trusted.
- **`tabperen` is `bool`** (Fortran `logical`), not float.
- **`Bad = -9999f` compared exactly** — no epsilon "fix".
- **Unported routines throw `NotImplementedException`** rather than returning zeros, so parity failures are loud.

---

## 4. Open blockers (in priority order)

1. **No golden corpus and no runnable oracle** — `bin/alexi_proc` is Linux ELF; on this host we need either a rebuilt Windows gfortran binary, a WSL/SSH run, or a checked-in captured corpus.
2. **No sample input files** — `INPUTS/ALEXI_INPUT/met_*.input` etc. do not exist in the repo, so the full CLI run path has never been exercised end-to-end.
3. **`Alexi.Golden.Tests` has only driver-helper tests** (GetDate/GetDgmt) — the parity harness (Fortran subprocess vs C# in-process) does not exist yet.
4. **Nothing is committed** — `.gitignore` now has `csharp/**/bin/`, `csharp/**/obj/`, `csharp/.vs/`; first commit still pending.
5. ~~InterpLAI compile errors~~ — **resolved**.

---

## 5. Next planned steps

**Sprint A — make the tree green (≈1 day)** — ✅ DONE except commit
1. ✅ Fixed `InterpLAI.cs`; ported all remaining modules (plan §5 items 5–13); full solution builds; 54 tests green.
2. ⏳ `.gitignore` extended — **commit still pending** (ask before committing).
3. ❌ CI workflow (build + `dotnet test`) — not added.

**Sprint B — build the oracle and the parity gate (Phase 0 completion)**
4. Produce a runnable Fortran oracle: rebuild `alexi_proc` with gfortran (`-ffixed-line-length-132 -g`) on a Linux host/WSL, or capture outputs once and check them in.
5. Build `tests/golden/`: 3–5 cases (clear-sky, cloudy, water/ice, non-converging, bad-input pixel) with `met_`/`sat_`/`nparm_`/`profile_`/`veg_us_*.input` inputs, expected stdout + `.LOG` + `*_YYYYDDD.dat` grids + per-cell numeric fingerprint.
6. Implement the parity harness in `Alexi.Golden.Tests` (Fortran as subprocess, C# in-process), tolerance `1e-3` absolute on fluxes.

**Sprint C — continue the port order (plan §5, modules 6–13)** — ✅ ALL DONE
7. ✅ `TwoSource` ← `ALEXI.f`.
8. ✅ `Water` ← `ALEXI_water.f`.
9. ✅ Clear-sky daily flux ← `USflux_clear.f`.
10. ✅ `Cover` completion ← `USflux_cover.f` + `landcover.f`.
11. ✅ `Alexi.IO` ← `USflux_run.f` I/O half.
12. ✅ `Input` ← `USflux_run.f` input half + `ListDirectedReader`.
13. ✅ `Driver` ← `USflux.f90` + `USflux_main.f` + CLI (4-arg contract + dir overrides + `--store-binary`).
Remaining CLI gap: pyalexi iteration flags (`--start-year/--end-doy`, `--point`, `--region`) not implemented.

**Sprint D — after parity is green only**
15. Config/data layout (plan §6): inject `AlexiPaths`, ship `store/landcover.txt` as an embedded resource with override.
16. Validation gate (plan §7): binary grid parity incl. `BAD` pattern and row flip; per-module debug dump (`ta1, ts1, tc1, ra2, rs2, rx2, z2, hn, psi`); property tests on `getfveg`/`psimhn`/`getsunzen`/`getxlwdn`; perf within 1.5× before parallelizing.
17. Post-parity wins (plan §8): `Parallel.For` over pixels, `Span<float>`/`ArrayPool` for the ~440 MB grids, structured logging, NuGet packaging (resolve the BSD/MIT/CC0 license conflict first).

---

## 6. How to verify current state

```powershell
cd csharp
dotnet build Alexi.sln        # 0 errors (legacy .sln; Alexi.slnx equivalent)
dotnet test Alexi.slnx        # 49 pass Alexi.Core.Tests + 5 pass Alexi.Golden.Tests
# note: bare `dotnet test` is ambiguous now that both .sln and .slnx exist
```

## 7. Caveats — do not trust yet

- Green unit tests ≠ parity. Nothing has been compared against Fortran output; all numeric assertions are self-derived.
- `FindAlbedoSoil`/`FindAlbedoVeg` convergence loops and `AlexiUtl.RunInit` are unvalidated against reference output.
- `DailyFluxCloudy` sets all daily fluxes to `BAD` — that mirrors the source's gap-fill behaviour, not a computed cloudy flux; confirm it is intended.
- The `lookup` alias decision (§3) is unconfirmed; any golden case touching `lookup_i/lookup_j` may encode undefined behaviour.
- `interp_LAI.f` is **not in the SCons build** — decide whether it belongs in the port at all before investing in it.
