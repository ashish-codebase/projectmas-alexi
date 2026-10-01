# ALEXI / projectmas_alexi — Master Handoff Document (v2)

**Repository:** `projectmas-alexi-master`
**Artifact:** an **unpacked conda package**, not a git repository
**Package:** `projectmas_alexi` v`0.1.0`, build `py37h47cba93_0`, `linux-64`, built 2020-08-31
**Upstream author of the Fortran:** Martha Anderson (USDA-ARS / UMD), 1996–2012
**Conda packaging author:** Mitchell Schull (NOAA), `mitch.schull@noaa.gov`

Every claim below was read out of the files in this tree. Where the existing `README.md` disagrees with the source, the discrepancy is listed in [§15 Corrections vs README.md](#15-corrections-vs-readmemd).

---

## Table of Contents

1. [What this tree actually is](#1-what-this-tree-actually-is)
2. [What ALEXI does (as implemented here)](#2-what-alexi-does-as-implemented-here)
3. [Directory map](#3-directory-map)
4. [Build chain](#4-build-chain)
5. [Runtime call graph](#5-runtime-call-graph)
6. [Input file specification](#6-input-file-specification)
7. [Output file specification](#7-output-file-specification)
8. [Include files and common blocks](#8-include-files-and-common-blocks)
9. [Per-file reference](#9-per-file-reference)
10. [Key equations as actually coded](#10-key-equations-as-actually-coded)
11. [Land cover classes — the real story](#11-land-cover-classes--the-real-story)
12. [Error codes](#12-error-codes)
13. [How to actually run it](#13-how-to-actually-run-it)
14. [Known bugs, landmines, and dead code](#14-known-bugs-landmines-and-dead-code)
15. [Corrections vs README.md](#15-corrections-vs-readmemd)
16. [Suggested next steps](#16-suggested-next-steps)

---

## 1. What this tree actually is

The top level is the **extracted contents of a conda `.tar.bz2` package**, not a source checkout:

```
projectmas-alexi-master/
├── README.md                 ← pre-existing reference doc (see §15 for its errors)
├── readme_v2.md              ← THIS FILE
├── bin/                      ← installed executables
│   ├── alexi_proc            ← compiled Fortran ELF, 341,712 B, linux-64
│   └── pyalexi               ← console-script shim (BROKEN, see §14.1)
├── lib/python3.7/site-packages/
│   ├── alexi-0.1-py3.7.egg   ← egg archive holding pyalexi.py
│   └── alexi.pth             ← 22-byte path hook into the egg
└── info/                     ← conda package metadata
    ├── about.json, files, git (empty), has_prefix, hash_input.json, index.json, paths.json
    ├── recipe/               ← the ONLY Fortran source in this tree
    └── test/run_test.py
```

`info/index.json` runtime dependencies:

```
python >=3.7,<3.8.0a0   numpy 1.17.*   libgcc-ng >=9.3.0   libgfortran-ng >=9,<10.0a0
```

`info/git` is a **zero-byte file** — no upstream git provenance recorded. `info/recipe/meta.yaml` records the original build host path `/PycharmProjects/projectMASalexi`, i.e. this package was built from a local working copy, not a pinned commit.

> **Consequence for handoff:** `bin/alexi_proc` is a prebuilt **Linux x86-64** binary. It will not run on Windows or macOS. To run ALEXI on this machine you must rebuild from `info/recipe/ALEXI_0.2/src/`.

Two source snapshots exist:

| Directory | Status |
|---|---|
| `info/recipe/ALEXI_0.1/src/` | Older snapshot. Same file names, no `landcover.f`, no free-form `USflux.f90` driver. **Not built.** |
| `info/recipe/ALEXI_0.2/src/` | **Active.** What `SConstruct.py` compiles and what produced `bin/alexi_proc`. |

All file references below are relative to `info/recipe/ALEXI_0.2/src/` unless stated otherwise.

---

## 2. What ALEXI does (as implemented here)

ALEXI (**A**tmospheric **L**and **EX**change **I**nverse) is a two-source (canopy + soil) energy-balance model that inverts for the **scaling sensible heat flux `Hn`**.

Given two geostationary (GOES) radiometric surface temperatures at two local times, plus surface meteorology and a boundary-layer potential-temperature profile, ALEXI:

1. Computes net radiation at both times from `SDN`, `XLWDN`, `albedo`, and estimated canopy/soil temperatures.
2. Splits the surface into a **canopy** source and a **soil** source using fractional cover `fc` derived from LAI.
3. Computes aerodynamic (`ra`), soil (`rs`) and leaf boundary-layer (`rx`) resistances from Monin–Obukhov similarity theory.
4. Advances the **convective boundary layer** from `t1` to `t2` using the Tennekes well-mixed slab model, and inverts the CBL heat budget for `Hn`.
5. Iterates with adaptive relaxation until `|ΔHn| < 0.1 W/m²`.
6. Scales the instantaneous flux at `t2` to a **daily total** using the daily-integrated shortwave (`sday`).

The inversion variable is defined in `ALEXI.f:findhn`:

```fortran
h1 = hn*t1/thrise      ! sensible heat flux at time 1  [W/m2]
h2 = hn*t2/thrise      ! sensible heat flux at time 2  [W/m2]
```

with `thrise = 3.0*3600 s` (Tennekes 1973 flux rise time). `Hn` is the flux the CBL would receive at `t = thrise`; the flux at any time grows linearly with time since the assumed sunrise offset.

**Domain:** `ilg = 1456`, `jlg = 625` (`USflux_parm.inc`) — 910,000 cells max, but the actual loop bound is the `npoints` command-line argument, not the grid size.

---

## 3. Directory map

```
info/recipe/
├── meta.yaml                 ← generated conda recipe (final, conda-build 3.19.2)
├── meta.yaml.template        ← hand-written recipe it was generated from
├── conda_build_config.yaml   ← python 3.7, numpy 1.17, gfortran, target linux-64
├── build.sh                  ← runs `$PYTHON setup.py install` (NOT scons directly)
├── setup.py                  ← actually invokes scons, then installs the python egg
├── SConstruct_testing.py     ← stale variant; builds a separate `landcover` binary
├── construct.yaml            ← conda "constructor" spec → name: pyalexi v1.0.0
├── LICENSE.txt
├── README.md                 ← leftover stub titled "# pyRTTOV" (see §14.9)
├── pyalexi/
│   ├── __init__.py           (empty)
│   └── pyalexi.py            ← the Python wrapper
├── ALEXI_0.1/src/            ← older snapshot, unused by the build
└── ALEXI_0.2/
    ├── src/                  ← ★ ACTIVE Fortran + SConstruct + Makefile
    ├── store/landcover.txt   ← class parameter table (7 usable rows)
    └── fluxes/               ← empty placeholder for binary output
```

`ALEXI_0.2/src/` also contains **committed build artifacts**: 15 `.o` object files and `.sconsign.dblite`. These are stale and should not be trusted as build state.

---

## 4. Build chain

The chain is longer than the old README implies:

```
conda-build
  └─ build.sh:  $PYTHON setup.py install
       └─ setup.py:
            cd <work>/ALEXI_0.1/src          # ← note: 0.1, not 0.2  (see §14.10)
            scons -Q --prefix=$PREFIX install
            scons -c
       └─ setuptools: installs pyalexi.py as an egg + console_scripts entry point
```

`ALEXI_0.2/src/SConstruct.py` is the real build definition:

```python
env.Append(F90FLAGS=['-g', '-fPIC', '-ffixed-line-length-132'])
alexi = env.Program(target='alexi_proc', source=[
    'USflux.f90','USflux_utl.f','USflux_rad.f','USflux_cover.f','USflux_run.f',
    'USflux_cloud.f','ALEXI.f','ALEXI_atmos.f','ALEXI_rad.f','ALEXI_utl.f',
    'ALEXI_water.f','pbl_read.f','sfc_read.f','landcover.f','USflux_main.f'])
env.Install(bin_path, [alexi])
```

Notes:

- **One executable only:** `alexi_proc`. `landcover` is *not* a separate program in the active build — `landcover.f` is linked in and called as a subroutine.
- `interp_LAI.f` is **not** in the SCons source list; it is only built by the legacy `Makefile` (`make all` → `alexi_proc landcover interp_LAI`).
- The `Makefile` is a legacy path (author MA, 03/23/2012). It expects an `./o` object directory and a `../bin` output directory that do not exist in this tree.
- SCons sets no `-mcmodel` / `-heap-arrays`; the Makefile link line uses `-mcmodel=medium`. See §14.7 for why memory matters.
- `SConstruct.py` depends on conda-build env vars `PREFIX`, `BUILD_PREFIX`, `SRC_DIR`, `FC`. It cannot run standalone without setting them.

**Rebuild on Linux/WSL:**

```bash
cd info/recipe/ALEXI_0.2/src
export FC=$(which gfortran) PREFIX=$PWD/_install BUILD_PREFIX=$PWD SRC_DIR=$PWD/..
scons -Q --prefix=$PREFIX install
# → _install/bin/alexi_proc
```

---

## 5. Runtime call graph

```
alexi_proc  <MDATE:YYYYDDD>  <tile:char*16>  <npoints:int>  <part:char*3>
│
└─ program main                        USflux_main.f
   ├─ getarg 1..4  → MDATE, t, npoints, part
   │
   ├─ call landcover(MDATE,t,npoints,part)                    landcover.f
   │    ├─ loadclasstable      → reads  <HOMEDIR>/store/landcover.txt
   │    ├─ READ   INPUTS/ALEXI_INPUT/veg_<MDATE>_<tile>_<part>.input      (unit 102)
   │    └─ WRITE  INPUTS/ALEXI_INPUT/veg_us<MDATE>_<tile>_<part>.input    (unit 103)
   │         per line: i j aleafv aleafn aleafl adeadv adeadn adeadl height xleaf z0eff dispeff rsmin
   │         (uses weighted_avg → canopyarch_lai)
   │
   └─ call USflux(MDATE,t,npoints,part)                      USflux.f90   ★ the pixel loop
        ├─ open 99  ALEXI_LOG/<MDATE>_<tile>_<part>.LOG
        ├─ call set_constants                                 ALEXI_utl.f
        ├─ open 102 INPUTS/ALEXI_INPUT/met_<MDATE>_...input
        ├─ open 110 INPUTS/ALEXI_INPUT/sat_<MDATE>_...input
        ├─ open 120 INPUTS/ALEXI_INPUT/nparm_<MDATE>_...input
        ├─ open 130 INPUTS/ALEXI_INPUT/profile_<MDATE>_...input
        ├─ open 400 INPUTS/ALEXI_INPUT/veg_us<MDATE>_...input
        ├─ open 103 OUTPUTS/ALEXI_OUTPUT/output_<MDATE>_...input
        │
        └─ do m = 1, npoints
             ├─ extract_input(m,ibad,MDATE)                  USflux_run.f
             │    ├─ reads one line from each of 102/110/120/130
             │    ├─ cover_props(ia,ja)                      USflux_cover.f  (reads unit 400)
             │    │    └─ getfveg / canopyarch
             │    └─ sets `clear` from iflag
             ├─ if clear:
             │    ├─ ALEXI(ia,ja,ierr,ibad,iter)             ALEXI.f        (land)
             │    │   or ALEXI_water(ia,ja,ierr,ibad,iter)   ALEXI_water.f  (inland water)
             │    │    ├─ runinit → checkinput, getsunzen, getpbltable, updatefc
             │    │    └─ loop: findhn → getnetrad, getsoilheat, findta, growPBL
             │    │         └─ getresistance / psimhn_ky     ALEXI_atmos.f
             │    └─ if converged: clear_day_proc            USflux_clear.f
             │         ├─ hourly_flux                        USflux_utl.f
             │         │    ├─ getsunzen, getradcomps        USflux_rad.f
             │         │    ├─ getradprops, getnetrad        ALEXI_rad.f
             │         │    ├─ getsoilheat                   ALEXI.f
             │         │    ├─ fao_PM → getnetrad_simple     USflux_rad.f
             │         │    └─ integrate  (×8)
             │         └─ daily_flux_clear_SDN
             │    else: cloudy_day_proc                      USflux_cloud.f
             │         ├─ hourly_flux
             │         └─ daily_flux_cloudy   (sets daily totals = BAD)
             └─ screen_output                                USflux_run.f
                  └─ WRITE unit 103: i j rnet2 xle2 h2 g2 xles2
```

`store_output_bin` (binary gridded output) exists in `USflux_run.f` but is **commented out at its only call site** (`USflux.f90` line 136). No binary output is produced by this build.

`pbl_read.f` and `sfc_read.f` are compiled and linked but **never called** from the active pipeline. They belong to the older gridded (CFSR/MENAE) driver. See §14.7.

---

## 6. Input file specification

All input files are **plain text, one line per pixel**, read sequentially by line number `m = 1..npoints`. Leading index columns are read but **not validated against `m`** — file order is what matters.

File names are built as `<prefix>_<MDATE>_<tile>_<part>.input` in `./INPUTS/ALEXI_INPUT/`.

### 6.1 `met_<MDATE>_<tile>_<part>.input` — surface meteorology (unit 102)

```
i1 i2 i3 i4  mvars(1) … mvars(10)
```

| # | Variable | Meaning | Units |
|---|---|---|---|
| — | i1,i2,i3,i4 | index columns (i3,i4 discarded) | — |
| 1 | `taobs1` | air temperature at t1 | **K** (converted to °C) |
| 2 | `taobs2` | air temperature at t2 | **K** |
| 3 | `ea1` | vapour pressure at t1 | **mbar** |
| 4 | `ea2` | vapour pressure at t2 | **mbar** |
| 5 | `pres1` | station pressure at t1 | **mbar** |
| 6 | `pres2` | station pressure at t2 | **mbar** |
| 7 | `xlwdn1` | downward longwave at t1 | W/m² |
| 8 | `xlwdn2` | downward longwave at t2 | W/m² |
| 9 | `wind1` | wind speed (measured at 30 m) | m/s |
| 10 | `wind2` | wind speed (measured at 30 m) | m/s |

Hard rules enforced in `extract_input`:
- `taobs1 > taobs2` → `badinput = .TRUE.` (cooling daytime air is rejected).
- `wind1 = wind2 = (mvars(9)+mvars(10))/2` — **the two wind columns are averaged; both times get the same value.**
- `xlwdn1/xlwdn2` read from the file are then **overwritten** by `call getxlwdn(ea,ta,1.0,xlwdn)` — the model recomputes downwelling longwave from humidity and ignores the input columns.
- Wind is rescaled from 30 m to the 50 m blending height, then **clamped to [3.0, 20.0] m/s** in `runinit`.

### 6.2 `sat_<MDATE>_<tile>_<part>.input` — satellite observations (unit 110)

```
i1 i2  svars(1) … svars(8)
```

| # | Variable | Meaning | Units |
|---|---|---|---|
| 1 | `trad1` | radiometric surface temperature at t1 | **K** |
| 2 | `trad2` | radiometric surface temperature at t2 | **K** |
| 3 | `sdn1` | net shortwave at t1 | W/m² |
| 4 | `sdn2` | net shortwave at t2 | W/m² |
| 5 | `r15` | observation time 1, **hours after sunrise** | h |
| 6 | `r55` | observation time 2, hours after sunrise | h |
| 7 | `theta` | satellite view zenith angle | rad |
| 8 | `xlai*10` | LAI × 10 (divided by 10 on read) | — |

Hard rules:
- `dthr = (trad2-trad1)/(r55-r15)`; `offset = (dthr*rsmin - dthr)*(r55-r15)`; `trad1 = trad1 - offset`. This is the **emissivity/resistance-corrected warming-rate adjustment** applied to the first radiometric temperature, scaled by `rsmin` from the `nparm` file.
- `trad1 > trad2` → `badinput = .TRUE.` (radiometric surface must warm between overpasses).
- `r55 - r15 <= 0` → `iflag = 2` (cloudy path).
- `xlai < 0` → `xlai = fc = height = z0 = disp = BAD`.
- `fc = 1 - exp(-0.5*xlai)`, floored at 0.01.
- `svars(9)` is referenced (`if (svars(9).eq.-9999.)`) but **only 8 values are read** → stale/uninitialised array element. See §14.3.

### 6.3 `nparm_<MDATE>_<tile>_<part>.input` — land/surface parameters (unit 120)

```
i1 i2  pvars(1) … pvars(6)
```

| # | Variable | Meaning | Used as |
|---|---|---|---|
| 1 | `rsoilv` | soil reflectance, visible | → `rsoilv` |
| 2 | `rsoiln` | soil reflectance, NIR | → `rsoiln` |
| 3 | `xlat` | pixel latitude | degrees |
| 4 | `xlong` | pixel longitude (0–360 convention handled by `getdgmt`) | degrees |
| 5 | `iclass` | land cover class | **IGNORED — hardcoded `iclass = 5`** |
| 6 | `rsmin` | minimum stomatal/soil resistance scalar | → `rsmin`, and drives the `dthr` correction |

> `iclass=5 !pvars(5)` at `USflux_run.f:97` is the single most important line in this codebase for interpreting outputs: **every land pixel is treated as class 5** regardless of the input file.

### 6.4 `profile_<MDATE>_<tile>_<part>.input` — boundary-layer sounding (unit 130)

```
i1 i2  rvars(1) … rvars(14)
```

`rvars(k)` = potential temperature (K) at the k-th level. Heights are **hardcoded** in `extract_input`:

```
0, 100, 300, 500, 700, 1000, 1400, 1800, 2200, 2600, 3000, 3500, 4000, 4500  [m]
```

These fill `zpbli(1:14)` / `thpbli(1:14)`. `ALEXI_parm.inc` declares `mli = 41` input levels; levels 15–41 are left undefined. `getpbltable` builds the lookup table used by `growPBL`.

### 6.5 `veg_<MDATE>_<tile>_<part>.input` — sub-pixel class mix (consumed by `landcover()`)

```
i  j  xntot  freq(1) … freq(nclass)  xlai*10
```

`nclass = 8`. `freq(k)` are fractional abundances of each class in the pixel. `xlai` is divided by 10; values `<0` or `>10` become `BAD`.

### 6.6 `veg_us<MDATE>_<tile>_<part>.input` — generated by `landcover()`, consumed by `cover_props()`

Written by `landcover()` (unit 103), read by `cover_props()` (unit 400):

```
i  j  aleafv aleafn aleafl adeadv adeadn adeadl height xleaf z0eff dispeff rsmin
format: 2i5, 8f13.5, 3f13.5
```

These **overwrite** whatever `cover_props()` computed from `landcover.txt`. The comment in `USflux_cover.f` is explicit: *"Currently, these data overwrite values computed from the table parameters."*

---

## 7. Output file specification

### 7.1 `OUTPUTS/ALEXI_OUTPUT/output_<MDATE>_<tile>_<part>.input` (unit 103)

Written by `screen_output` for every pixel that reaches it, format `2I6,5f10.1`:

```
ai  aj  rnet2  xle2  h2  g2  xles2
```

| Column | Meaning | Units |
|---|---|---|
| `ai`,`aj` | index pair carried from the input line (not grid i,j) | — |
| `rnet2` | net radiation at t2 | W/m² |
| `xle2` | **total** latent heat flux at t2 | W/m² |
| `h2` | **total** sensible heat flux at t2 | W/m² |
| `g2` | soil heat flux at t2 | W/m² |
| `xles2` | **soil** latent heat flux at t2 | W/m² |

Canopy latent `xlec2 = xle2 - xles2` must be derived by the user.

### 7.2 `ALEXI_LOG/<MDATE>_<tile>_<part>.LOG` (unit 99)

Run summary at the end of `USflux()`:

```
Bad input: nbad | Water or ice: nwater | Cloudy pixel: ncloud
Converged: nconv | Did not converge: nfail | Total points: ntot
```

Plus `write(99,1003)` diagnostic lines (`TC1 TS1 TC2 TS2 HC2 HS2 ZEN1 ZEN2 RNET1 TA1 TA2 TAOBS1 TAOBS2 TH1 TH2`) and the daily totals block (`RNDAY HDAY GDAY EDAY ESDAY ECDAY EREFDAY SDAY`). The daily-totals block sits under a commented-out `elseif` — see §14.4.

### 7.3 Screen (unit 6)

Verbose per-pixel diagnostics gated on `writeme`, which `extract_input` sets `.FALSE.` for every pixel. In practice the screen shows only the banner and stray debug writes.

### 7.4 Daily-integrated products (computed, mostly not written)

`daily_flux_clear_SDN` computes `rnday, rnsday, rncday, gday, eday, esday, ecday, hday, hsday, hcday, sday, aparday, acday` plus `fpet`, `fsdn`. Integrated fluxes are in **MJ/m²/day** (`integrate()` does `Σf / xn * 3600 * 24 / 1e6`).

```fortran
fsdn  = sdn2 / sday          ! instantaneous-to-daily shortwave ratio
eday  = xle2 / fsdn          ! daily ET  [MJ/m2/day]
gday  = 0.0                  ! G assumed to integrate to zero over 24 h
hday  = rnday - eday - gday
esday = (xles2/xle2) * eday
hsday = rnsday - gday - esday
ecday = eday - esday
hcday = rncday - ecday
fpet  = xle2 / eref2         ! FAO reference-ET ratio, floored at 0
```

---

## 8. Include files and common blocks

| File | What it actually contains |
|---|---|
| `ALEXI_parm.inc` | `mli=41` (input PBL levels), `ml=8000` (interpolated levels), `mt=8000` (lookup table size), `BAD=-9999.0`, `SUCCESS=0`. **Nothing else** — no κ, no cp, no grid dims. |
| `ALEXI.inc` | Fixed-form declaration of the 27 ALEXI common blocks + their `logical` type statements. |
| `ALEXI_f90.inc` | Free-form (`&` continuation) duplicate of `ALEXI.inc`, used by `USflux.f90`. The two must be kept in sync by hand. |
| `USflux_parm.inc` | `kz=30, ilg=1456, jlg=625` (MENAE domain); `kx=1440, ky=600, kt=8` (GETD input grid); `sx=1440, sy=720` (clear-sky INSOL grid); `nohrc=24, nohr=24`; `nclass=8`, `ICE=WATER=RIPARIAN=1`; `nclasm=6`, `nbin=20`; NCAR-graphics sizes. |
| `USflux_grids.inc` | 2-D `real*4` arrays in `/canopy1/`, `/canopy2/`, `/canopy3/`. `CLAT, CLON, DX, DY, MINLAT…` are **uninitialised variables**, never set anywhere in this tree. |
| `USflux_dir.inc` | `HOMEDIR="./ALEXI/"`, `INDIR="./ALEXI/inputs/"`, `TABDIR="./ALEXI/store/"`, `OUTDIR="./ALEXI/fluxes/"`, `INPUT="./ALEXI/inputs/"` — all relative to CWD. |
| `USflux.inc` | US-layer common blocks: `/cover/`, `/dayflux*/`, `/hrdata*/`, `/lookup/`, `/stress/`, `/sun/`, `/usflags/`. Declares `CHARACTER(LEN=50) tabdesc`. |
| `date.inc` | `/canopy0/ JDATE, MDATE, CDATE, YDATE` (`integer JDATE, MDATE`; `character*8 CDATE`; `character*17 YDATE`). |
| `flag.inc` | `integer*8 filg=1456, fjlg=625`; `common/fflag/flag(filg,fjlg)`. |

**Directory inconsistency:** `USflux_dir.inc` declares `./ALEXI/inputs/` and `./ALEXI/fluxes/`, but `USflux.f90` and `landcover.f` hardcode `./INPUTS/ALEXI_INPUT/`, `./OUTPUTS/ALEXI_OUTPUT/`, `./ALEXI_LOG/`. Only `TABDIR`/`HOMEDIR` are honoured (`load_tables` reads `TABDIR//'landcover.txt'`; `loadclasstable` reads `HOMEDIR//'store/landcover.txt'`). See §13.

---

## 9. Per-file reference

All paths relative to `info/recipe/ALEXI_0.2/src/`.

### 9.1 Physics core

#### `ALEXI.f` — the inverse solver

| Subroutine | Purpose |
|---|---|
| `alexi(ia,ja,ierr,ibad,iter)` | Entry point. `delhmx = 0.1 W/m²`. Outer loop over `fc` (currently disabled). Inner `do while(.not.converged)`: `findhn` → relax `hn = hn + f*delthn` with `f = 0.25/(1+int(iter/20))`, `iter < 300`. On success caches `hn0`, `psi0` as the first guess for the **next pixel**. |
| `findhn(hn,hnnew,converged,stopiter,ierr)` | One forward-model evaluation: `getnetrad` ×2 → `getsoilheat` ×2 → `findta` at t1 → `findta` at t2 → `growPBL` → returns `hnnew`. Caps `h2new ≤ rnet2 − g2`. Sets `ierr = 30` if `hnnew ≤ 0`. |
| `growPBL(hnnew,stopiter,ierr)` | CBL growth inversion (Tennekes slab): integrates `H` over `[t1,t2]` against the potential-temperature profile to obtain the required `Hn`. |
| `growPBL2(...)` | Alternate/legacy CBL integrator, **not called** from `findhn`. |
| `findta(trad,ta,ts,tc,rn,rnsoil,rndiv,h,hs,hc,xle,xles,xlec,g,ra,rs,rx,wind,rhocp,pres,ea,tac,stopiter,ierr,itype)` | Solves the canopy/soil temperature split consistent with `trad` and `h`. |
| `findta_parallel(...)` | Variant, **not called**. |
| `getsoilheat(rnsoil,g,tloc)` | `G/Rsoil` parameterisation as a function of time since sunrise. |
| `getsoilheat_old(rnsoil,g)` | Legacy, **not called**. |
| `checksoln(rnsoil,g,ts,trad,converged,ierr)` | Solution plausibility test; **call is commented out** in `alexi()`. |

State carried across pixels: `hn0` and `psi0` are `save`d and updated on every converged pixel, so each pixel's starting guess is the previous pixel's solution. Results are therefore **order-dependent**.

#### `ALEXI_atmos.f` — resistances and canopy architecture

| Subroutine | Purpose |
|---|---|
| `getresistance(h,wind,ta,ts,tc,ra,rs,rx,rhocp)` | Monin–Obukhov, two-pass `u*` (buoyancy first, then stability-corrected). Water/bare-soil branch uses `z0h = z0*exp(−κ(4·Re*^0.15 − 5))` and sets `rs = 0`, `rx = 10000`. Vegetated branch: `ra = (xlog1−ψm)(xlog1−ψh)/(0.16·wind)`, `rs = 1/(0.0025·ΔT^0.33 + 0.012·us)` with `ΔT = ts − ta`, `rx = 180·sqrt(xl/udz)/xlai`. |
| `psimhn_ky(zdla,psima,psih)` | Kader & Yaglom (1990) stable/unstable ψ_m, ψ_h. **Active.** |
| `psimhn(zdla,psima,psih)` | Classical Dyer/Hogström ψ. **Not called.** |
| `canopyarch` | Sets `z0`, `disp` from `fc` and `height` for the current pixel. |

Note: the NDVI-dependent variants of `val1`/`val2` in the `rs` formulation are present but commented out; the constants are fixed at `0.0025` and `0.012`.

#### `ALEXI_rad.f` — radiation

| Subroutine | Purpose |
|---|---|
| `getsunzen(xlat,xlong,stdlng,doy,year,ftime,zen)` | Solar zenith angle. |
| `getradprops(zen,albedo,taubtv,taubtn,clumps)` | Canopy albedo and transmittance (Goudriaan 1988 two-stream, VIS/NIR). |
| `getnetrad(ts,tc,xlwdn,sdn,albedo,taubtv,taubtn,rnet,rnsoil,rndiv,swup,xlwup)` | Net radiation **above canopy** and **above soil**. |
| `getnetradnight(...)` | Nighttime variant. |
| `getxlwdn(ea,ta,fclear,xlwdn)` | Downwelling longwave from vapour pressure (Brutsaert-type). |
| `find_albedo_soil(tloc2)` / `find_albedo_veg(tloc2)` | Albedo estimation helpers; present but not on the main path (`albedo2` is set to `BAD` in `extract_input`). |

#### `ALEXI_water.f` — inland-water solver

| Subroutine | Purpose |
|---|---|
| `alexi_water(ia,ja,ierr,ibad,iter)` | Same relaxation loop as `alexi()`, no canopy/soil split. |
| `findflux_water(hn,hnnew,converged,stopiter,ierr)` | Water version of `findhn`. |
| `findta_water(...)` | `ta = trad − h·ra/rhocp`; `xle = rn − h − g`. |
| `getsoilheat_water(rnet,g)` | **`g = 0.55*rnet`** (water heat storage). |
| `getnetrad_water(trad,xlwdn,sdn,rnet,swup,lwup)` | **`albedo = 0.1`, `emissivity = 0.99`** fixed. |

Selected by `USflux.f90` when `iswater_inland` is true — but `extract_input` sets `iswater_inland = .FALSE.` unconditionally, so in this build **`ALEXI_water` is unreachable** unless that logic changes.

#### `ALEXI_utl.f` — constants, init, validation

| Subroutine | Purpose |
|---|---|
| `set_constants` | `xk = 0.4`, `pi`, `thrise = 3*3600 s`, `cp = 1010 J/kg-K`, `z1 = 50 m`, `fcbare = 0.00`, clear-sky partition defaults, `hn0 = 80`, `psi0 = 0`. |
| `getprofile(z1)` | Interpolate the input potential-temperature profile to 1 m levels. **Call commented out in `runinit`.** |
| `getpbltable(z1,ierr)` | Build `tabtheta`/`tabz2` lookup tables used by `growPBL`. |
| `runinit(badinput,ibad,writeme)` | Per-pixel init: `checkinput`, clear-sky partitioning, `getsunzen` at t1/t2, sunrise offset `dtloc` (**hardcoded to 1.50 h**), emissivity coefficients (`eleaf = 0.97`), `updatefc`, wind rescale 30 m → 50 m plus clamp [3, 20] m/s, air density, `getpbltable`. |
| `updatefc` | Recomputes `esfc = aem·fc² + bem·fc + emsoil`. |
| `checkinput(flag,ibad,w)` | Range checks on ALEXI-level inputs. |
| `checkvalue(valname,itag,value,xmin,xmax,badinput,w,ibad)` | Generic range validator; assigns `ibad` codes 22–49. |
| `ALEXI_errorcode(writeme,ierr)` | Prints the convergence message. |

### 9.2 US driver layer

#### `USflux_main.f` — `program main`
Parses 4 command-line arguments, prints the banner, calls `landcover()` then `USflux()`. That is the whole program.

#### `USflux.f90` — `subroutine USflux(MDATE,t,npoints,part)`
**This is the pixel-loop driver**, not an alternative physics solver (see §15). Opens all input/output units, loops `m = 1..npoints`, dispatches clear/cloudy and land/water, tallies `ntot/nbad/nwater/ncloud/nconv/nfail`.

Dead branch: after a non-converged pixel, `go to 1000` executes before `call cloudy_day_proc(ibad)`, so the fallback cloud path for failed pixels is unreachable.

#### `USflux_run.f` — I/O (1181 lines, the largest file)

| Subroutine | Status |
|---|---|
| `extract_input(m,ibad,MDATE)` | **Active.** Reads the four input lines, builds `tloc1/tloc2` from `r15/r55 + dgmt`, applies the `dthr`/`rsmin` correction to `trad1`, sets flags. |
| `getdgmt(xnlon,dgmt,stdlng)` | **Active.** Standard longitude and local-time offset from a 0–360 longitude. |
| `screen_output(ia,ja)` | **Active.** Writes unit 103 and the diagnostics. |
| `store_output_bin(i,j,ierr,ibad,iflag,iter)` | Compiled, **call commented out**. |
| `open_output`, `binopen`, `binopen2`, `binwrite`, `binread`, `binread_swap`, `byteswapr4` | Binary-grid plumbing for `store_output_bin`; dormant. |
| `writepoint(ia,ja)` | Text-point writer, not called on the active path. |
| `extract_hourly_input_nldas(ia,ja,dgmt)` | NLDAS hourly time-series reader; **call commented out** in `extract_input`. |

#### `USflux_utl.f` — hourly energy balance and daily integration

| Subroutine | Purpose |
|---|---|
| `hourly_flux` | 24-hour loop: `getsunzen` → `checkvalue` → `getradcomps` → `getradprops` → proxy `ts`/`tc` from quadratic temperature departures (`getTdepart2`) → `getnetrad` → `getsoilheat` → `fao_PM`; then 8 × `integrate` to daily totals. |
| `fao_PM(eref,tloc,sdn,zen,ta,pres,ea,wind)` | FAO-56 Penman–Monteith grass reference ET, returns **W/m²**. `albedo = 0.23` hardcoded. |
| `integrate(f,fint,tstrt,tend,t,nhr)` | Sums hourly values and converts to MJ/day. |
| `getTdepart2(xm,xb,xc,trad,ta,t2,tr,ts)` | Fits the quadratic `ΔT(t)` departure used to build hourly `ts`/`tc` proxies. |
| `checkUSinput(flag,ibad,w)` | US-layer input validation; **call commented out** in `extract_input`. |
| `setPBLheights` | Builds a 41-level, 200 m-spaced `zpbli`. **Not called** — `extract_input` hardcodes 14 levels instead. |

Two in-source comments flag known approximations: `! CHANGE: clump should be clumps(ihr) !!!` (twice in `hourly_flux`) and `! CHANGE THESE DIMENSIONS!` in `integrate`.

#### `USflux_clear.f` — clear-sky daily aggregation

| Subroutine | Purpose |
|---|---|
| `clear_day_proc(ibad)` | `hourly_flux` then `daily_flux_clear_SDN`. |
| `daily_flux_clear_SDN` | **Active.** Scales t2 instantaneous fluxes to daily totals by `fsdn = sdn2/sday`. Sets `gday = 0`. Computes `fpet = xle2/eref2`. Sets `fawsfc = BAD`, `fawrz = BAD` — soil/canopy water-stress factors are **not supported** on this path. |
| `daily_flux_clear_EF` | Labelled **OLD — UNSUPPORTED** in its own header. Uses evaporative-fraction scaling and Jury & Tanner potential soil evaporation. Not called. |

#### `USflux_cloud.f` — cloudy path

| Subroutine | Purpose |
|---|---|
| `cloudy_day_proc(ibad)` | `hourly_flux` then `daily_flux_cloudy`. |
| `daily_flux_cloudy` | Sets `eday, hday, ecday, esday, hcday, hsday = BAD`. Header: *"Do not compute daily fluxes for cloudy pixels - gap-fill in post-processing."* |

#### `USflux_cover.f` — canopy properties

| Subroutine | Purpose |
|---|---|
| `load_tables(MDATE)` | Reads `TABDIR//'landcover.txt'` into the `/lookup/` tables; then hardcodes 14 `tabdesc` names. |
| `cover_props` | Derives `xlai`, `fg`, `clump*`, `height`, absorptivities from the tables — then **overwrites** `aleafv…adeadl, height, xl, z0, disp, rsmin2` by reading one line from unit 400 (`veg_us…input`). |
| `getfveg(clump0,xlai,fveg)` | Solves `exp(−0.5·clump0·xlai) = fv·exp(−0.5·xlai/fv) + (1−fv)` by brute-force scan over `fv = 0.01…1.00`; clamps `fveg` to [0.1, 1.0]. |

`perennial` is read from the table then **forced false**: `perennial=.FALSE.  ! THIS FUNCTIONALITY IS NOT ROBUSTLY IMPLEMENTED YET`.

#### `USflux_rad.f` — US-layer radiation utilities

| Subroutine | Purpose |
|---|---|
| `getradcomps(sdn,zen,fclear)` | Splits `SDN` into VIS/NIR and direct/diffuse following Cupid's `RADIN4`. Clearness `fclear = sdn/(potvis+potnir)`. Nighttime branch zeroes `sdn`. |
| `getnetrad_simple(sdn,ta,fclear,zen,albedo,rnet)` | Net radiation from `SDN` + air temperature only, used inside `fao_PM`. Sky emissivity: clear = Swinbank `9.2e-6·T⁴`; cloudy = Monteith & Unsworth blend on `fclear`. |

#### `landcover.f` — sub-pixel class averaging

| Subroutine | Purpose |
|---|---|
| `landcover(MDATE,tile,npoints,part)` | Reads `veg_…input`, writes `veg_us…input`. |
| `weighted_avg(...)` | Abundance-weighted mean of height, leaf size, absorptivities, `rsmin`; combines roughness via `z0eff = (50 − dispeff)/exp(1/sqrt(zsum))`. Water class uses `z0w = 0.00035 m`. |
| `canopyarch_lai(fc,hc,z0,disp)` | Simplified Massman form with hardcoded `dispdh = 2/3`, `z0dh = 1/8`; bare soil `z0s = 0.005`. The full Massman–Weil expressions are present as comments. |
| `loadclasstable` | Reads `HOMEDIR//'store/landcover.txt'` into `/lookup/` (`alv…adl, hmin, hmax, xl, rs`). |
| `getdate(MDATE,year,doy)` | `doy = mod(MDATE,1000)`, `year = MDATE/1000 − doy/1000`. |
| `laiopen`, `lairead`, `laiswapr4` | Legacy binary LAI grid readers pointing at the hardcoded absolute path `/data/data123/chain/4KM/GBIM/inputs/MLAI`. `laiopen` is commented out at its only call site. |

### 9.3 I/O and readers

**`pbl_read.f`** — `pbl_read(MDATE)` reads CFSR upper-air profiles and converts potential temperature to actual temperature over the `kz = 30` level MENAE domain; `fill_hgt_domain(above_hgt)` fills missing heights. **Compiled but never called** — the active path takes the sounding from the `profile_…input` text file.

**`sfc_read.f`** — `sfc_read(MDATE)` reads gridded CFSR/WRF surface fields into `/canopy1/` and sets up `navlat`/`navlon`; `arraynav` handles navigation. **Compiled but never called.**

**`interp_LAI.f`** — standalone `PROGRAM INTERP_LAI`. Reads `INPUT//'interp.in'` containing `file1, file2, outfile, f, ilg, jlg`, then linearly interpolates two unformatted 4-byte direct-access grids: `val = v1 + f*(v2 − v1)`, clamped at 0, `BAD` if either input `< −0.1`. Guarded by `ilgmx = jlgmx = 2000`. Built only by the legacy `Makefile`.

### 9.4 Standalone tools

| Tool | Build path | Status |
|---|---|---|
| `alexi_proc` | `SConstruct.py` | Built, shipped in `bin/` |
| `interp_LAI` | `Makefile` only | Not built by conda |
| `landcover` (standalone) | `Makefile` / `SConstruct_testing.py` | Not built by conda; now a subroutine inside `alexi_proc` |

### 9.5 Python layer

`info/recipe/pyalexi/pyalexi.py`:

```python
def alexi(year, doy, grid_name, npoints, part):
    date = "{}{:03d}".format(year, doy)
    subprocess.call(['alexi_proc', date, grid_name, npoints, part])
```

`main()` parses `--year`, `--doy`, `--grid_name` (`nargs='*'`, default `None`), `-c/--corr`, `-d/--debug`, then calls `alexi(args.year, args.doy, args.grid_name, npoints=None, part=0)` exactly once.

`bin/pyalexi` (the installed console script) is generated from the entry point `pyalexi=pyalexi.pyalexi:alexi` and ends with `sys.exit(alexi())`.

---

## 10. Key equations as actually coded

**Energy balance (instantaneous, at t2):**
```
Rn     = LE + H + G
Rnsoil = LEs + Hs + G
Rnc    = LEc + Hc
```

**Inverse variable:** `H(t) = Hn · t / thrise`, with `thrise = 10800 s`.

**Adaptive relaxation (`ALEXI.f`):**
```
ΔHn = Hn_new − Hn
if |ΔHn| > 0.1 W/m² and iter < 300:
    f  = 0.25 / (1 + int(iter/20))        ! 0.25 → 0.0125 across 300 iterations
    Hn = Hn + f·ΔHn
else if |ΔHn| ≤ 0.1: converged
```

**Resistances (`ALEXI_atmos.f`, vegetated):**
```
Ra = (ln((z−d)/z0) − ψm)(ln((z−d)/z0) − ψh) / (0.16·wind)
Rs = 1 / (0.0025·(Ts−Ta)^0.33 + 0.012·u_s)     if Ts−Ta > 1, else 1/(0.0025 + 0.012·u_s)
Rx = 180·sqrt(xl/u_dz) / LAI
```
`z = refhtw = zta = 50 m` (blending height); `xl` = leaf size (m). Note `ΔT` uses **`Ts − Ta`**, not `Ts − Tc` — the source comments: `! TC gets crazy at low fc`.

**Bare soil / water branch:**
```
z0h = z0·exp(−κ(4·Re*^0.15 − 5))
Ra  = (ln((z−d)/z0h) − ψh)/(u*·κ) ;  Rs = 0 ;  Rx = 10000
```

**Fractional cover from LAI:** `fc = 1 − exp(−0.5·LAI)`

**Canopy gap fraction (`getfveg`):** solve `exp(−0.5·clump0·LAI) = fv·exp(−0.5·LAI/fv) + (1 − fv)` for `fv`.

**Canopy architecture (`canopyarch_lai`):** `d = (2/3)·h`, `z0 = (1/8)·h`, `z0 ≥ 0.005 m`.

**FAO-56 Penman–Monteith (`fao_PM`, hourly, returns W/m²):**
```
λ    = 2.501 − 0.00237·Ta                       [MJ/kg]
γ    = 0.1·P·cp/(λ·0.622)                       [kPa/°C]   (P in mbar → kPa via ×0.1)
es   = 0.6108·exp(17.2694·Ta/(237.3+Ta))        [kPa]
s    = 17.2694·237.3·es/(237.3+Ta)²
D    = es − 0.1·ea                              (ea in mbar → kPa)
G    = 0.1·Rn
E_mm = [0.408·s·(Rn−G) + γ·37·u·D/(Ta+273.15)] / [s + γ·(1 + 0.24·u)]
E_W  = E_mm·λ·1e6/3600
```
⚠️ The aerodynamic terms use **`37·u` and `0.24·u`**, not the FAO-56 daily standard `900·u` and `0.34·u`. This corresponds to the hourly formulation, but verify against your reference before publishing ET ratios. Reference-surface albedo is hardcoded `0.23`.

**Daily scaling (`daily_flux_clear_SDN`):** `E_day = LE(t2)·sday/SDN(t2)`, `G_day = 0`.

**Water solver:** `albedo = 0.1`, `emissivity = 0.99`, `G = 0.55·Rn`.

---

## 11. Land cover classes — the real story

There are **three inconsistent class systems** in this tree. This is the most confusing part of the codebase.

### 11.1 The parameter table: `ALEXI_0.2/store/landcover.txt`

`nclass = 8`. The file has a header row plus 8 data rows — but the 8th row is malformed (it sits after the table with different whitespace and is labelled class **5**, duplicating Cropland as "Wetland"). Because `loadclasstable` reads exactly `nclass = 8` rows and `load_tables` reads until EOF keyed on `iclass`, that trailing row **overwrites class 5's parameters with the Wetland values**.

| Class | aleafv | aleafn | aleafl | adeadv | adeadn | adeadl | beta (LUE) | hmin (m) | hmax (m) | xl (m) | MODIS | Peren | fcmin | Description |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | −9999 | −9999 | −9999 | −9999 | −9999 | −9999 | 100 | −9999 | −9999 | −9999 | 0 | 0 | −9999 | Water |
| 2 | 0.86 | 0.37 | 0.95 | 0.84 | 0.61 | 0.95 | 150 | 5.0 | 10.0 | 0.10 | 0 | 0 | −9999 | Forest |
| 3 | 0.85 | 0.37 | 0.95 | 0.72 | 0.44 | 0.95 | 40 | 1.0 | 1.0 | 0.02 | 0 | 0 | −9999 | Shrubland |
| 4 | 0.82 | 0.28 | 0.95 | 0.42 | 0.04 | 0.95 | 40 | 0.1 | 2.0 | 0.02 | 0 | 0 | −9999 | Grassland |
| 5 | 0.83 | 0.35 | 0.95 | 0.49 | 0.13 | 0.95 | 40 | 0.1 | 2.0 | 0.05 | 0 | 0 | −9999 | Cropland |
| 6 | 0.82 | 0.57 | 0.95 | 0.92 | 0.80 | 0.95 | 40 | 1.0 | 1.0 | 0.02 | 0 | 0 | −9999 | Bare_Ground |
| 7 | 0.84 | 0.37 | 0.95 | 0.58 | 0.26 | 0.95 | 100 | 6.0 | 6.0 | 0.02 | 0 | 0 | −9999 | Urban_and_Built-Up |
| 5 ⚠ | 0.83 | 0.35 | 0.95 | 0.49 | 0.13 | 0.95 | 40 | 0.1 | **0.2** | 0.05 | 0 | 0 | −9999 | Wetland (malformed duplicate row) |

Class 8 is therefore **never populated** — `tabaleaf(8,*)` etc. remain uninitialised.

### 11.2 Hardcoded descriptions in `load_tables`

`USflux_cover.f` assigns `tabdesc(1)`…`tabdesc(14)` with MODIS-style names:

```
1  Water (and Goodes Interrupted Space)     8  Wooded Grassland
2  Evergreen Needleleaf Forest              9  Closed Shrubland
3  Evergreen Broadleaf Forest              10  Open Shrubland
4  Deciduous Needleleaf Forest             11  Grassland
5  Deciduous Broadleaf Forest              12  Cropland
6  Mixed Cover                             13  Bare Ground
7  Woodland                                14  Urban and Built-Up
```

But `tabdesc` is dimensioned `tabdesc(nclass)` = `(8)`. **Writing indices 9–14 is an out-of-bounds write** past the end of the `/lookup/` common block (§14.2).

### 11.3 What actually reaches the physics

`extract_input` sets `iclass = 5` unconditionally. So:

- `load_tables` / `cover_props` table lookups resolve to class 5 — which, after the malformed overwrite, are the **Wetland** values (`hmax = 0.2 m`).
- …and `cover_props` then **overwrites** `aleafv…adeadl, height, xl, z0, disp, rsmin2` from the per-pixel `veg_us…input` file produced by `landcover()`.
- `landcover()` itself uses `loadclasstable`'s copy of the same table, weighted by the `freq(1..8)` abundances in `veg_…input`.

**Net effect:** the class table matters mainly through `landcover()`'s weighted averaging; the `iclass` column in the `nparm` file is decorative.

---

## 12. Error codes

### 12.1 `ierr` — ALEXI solver status

Defined implicitly across `ALEXI_parm.inc`, `ALEXI.f`, `ALEXI_water.f`, and printed by `ALEXI_utl.f:ALEXI_errorcode`.

| Code | Meaning | Where set | Message printed |
|---|---|---|---|
| 0 | **CONVERGED** | `alexi()` | `*** CONVERGED ***` |
| 1 | `Hn` did not converge within 300 iterations | `alexi()` | `*** HN did not converge. Bail ***` |
| 2 | Bad input at `runinit` | `alexi()` | `*** Bad input. Bail ***` |
| 3 | Exceeded max PBL profile layer | `growPBL` | `*** Exceeded max PBL profile layer. Bail ***` |
| 9 | `Rnet2 < 0` | `findhn` | *(commented out)* |
| 11 | `Ts < Tc` for all `fc` | *(commented out)* | *(commented out)* |
| 12 | `G/Rnsoil < 0.1` | `checksoln` | *(commented out)* |
| 13 | `Ts − Trad > 15` | `checksoln` | *(commented out)* |
| 30 | `Hn_new ≤ 0` | `findhn` | *(no message in `ALEXI_errorcode`)* |

Codes 4–8 and 10 exist in the switch but print nothing.

### 12.2 `ibad` — input-validation code (assigned by `checkvalue` calls)

Observed assignments: `22` XLAI, `37` PRES/EA, `38` TA, `39` WIND, `40` SDN, `41` XLWDN, `42` ICLASS, `44` FG, `47` all insolation slots bad, `49` hourly-loop failure in `hourly_flux`.

### 12.3 `iflag` — per-pixel status (`USflux_run.f` comment block)

```text
0  good input data
1  bad input data (not cloud)
2  cloudy conditions
3  water / ice
4  ALEXI attempted but failed to converge   (assignment commented out in USflux.f90)
```

---

## 13. How to actually run it

`bin/alexi_proc` is Linux-only and the code resolves **relative paths against the current working directory**, so running it requires this scaffold:

```text
workdir/                          ← cd here before running
├── ALEXI/
│   └── store/
│       └── landcover.txt         ← required (HOMEDIR//'store/landcover.txt' and TABDIR//'landcover.txt')
├── INPUTS/
│   └── ALEXI_INPUT/
│       ├── met_<YYYYDDD>_<tile>_<part>.input
│       ├── sat_<YYYYDDD>_<tile>_<part>.input
│       ├── nparm_<YYYYDDD>_<tile>_<part>.input
│       ├── profile_<YYYYDDD>_<tile>_<part>.input
│       └── veg_<YYYYDDD>_<tile>_<part>.input
├── OUTPUTS/
│   └── ALEXI_OUTPUT/             ← must already exist; open() will not create it
└── ALEXI_LOG/                    ← must already exist
```

Copy `info/recipe/ALEXI_0.2/store/landcover.txt` into `workdir/ALEXI/store/` (and fix the malformed row first — §14.6).

**Direct invocation:**

```bash
cd workdir
../bin/alexi_proc 2017187 us 1000 0
#                 │       │   │     └ part: character*3, appended to every filename
#                 │       │   └─────── npoints: integer, number of input LINES to process
#                 │       └─────────── tile label: character*16, free text
#                 └─────────────────── MDATE = YYYYDDD, read with read(arg1,'(i7)')
```

`npoints` is a **line count**, not a grid index. It must not exceed the number of lines in each input file, and the four input files must have matching line counts and matching pixel ordering.

**Via the Python wrapper** (only works after fixing §14.1):

```bash
pyalexi --year 2017 --doy 187 --grid_name us
```

`--grid_name` is `nargs='*'` and is interpolated with `'{}'.format(...)`, so a multi-element list becomes a Python list literal as a shell argument. Use exactly one token.

---

## 14. Known bugs, landmines, and dead code

Ordered by how likely they are to bite you.

### 14.1 `pyalexi` CLI is broken as shipped
`bin/pyalexi` ends in `sys.exit(alexi())`, but `alexi()` is `alexi(year, doy, grid_name, npoints, part)` — five required positional arguments. Running `pyalexi` raises `TypeError`. The entry point should target `pyalexi.pyalexi:main`. Additionally `main()` calls `alexi(..., npoints=None, part=0)`, so `subprocess` receives the literal string `'None'` for `npoints`, and `read(arg2,'(i10)') npoints` fails. `num_parts`, `parts`, `doys`, `years`, `--grid_name` and `--corr` are computed but never used.

### 14.2 Out-of-bounds write in `load_tables`
`tabdesc` is dimensioned `(nclass)` = `(8)`, but `USflux_cover.f` assigns `tabdesc(9)`…`tabdesc(14)`. This writes past the end of the `/lookup/` common block. It happens to be the last item in the block, so it usually goes unnoticed, but it is undefined behaviour and will corrupt memory under any compiler that sizes the block exactly.

### 14.3 Uninitialised array read in `extract_input`
`svars` is declared `real svars(8)` and only 8 values are read, yet the code tests `if (svars(9).eq.-9999.)`. `fg` is then set to `1.0` regardless. Reading element 9 is undefined.

### 14.4 Dead branch in `screen_output`
The daily-totals block (`RNDAY HDAY GDAY EDAY ESDAY ECDAY EREFDAY SDAY`) sits under a commented-out `elseif (writeme.and..not.clear)`; it is now part of the main `if (writeme.and.clear.and.converged)` block. Combined with `writeme` being `.FALSE.` for every pixel, **no per-pixel diagnostics reach either the screen or the LOG file** in normal operation.

### 14.5 `cover_props` argument mismatch
`extract_input` calls `call cover_props(ia,ja)` but `USflux_cover.f` defines `subroutine cover_props` with **no arguments**. It compiles only because the call site is fixed-form Fortran with no explicit interface; the extra arguments are silently ignored.

### 14.6 `landcover.txt` malformed trailing row
The "Wetland" row is labelled class 5 and overwrites Cropland's parameters (`hmax` 2.0 → 0.2 m). Class 8 is never populated. Fix the file before trusting any canopy-height output.

### 14.7 `pbl_read.f` and `sfc_read.f` are dead weight
Compiled and linked, never called. They reference the gridded MENAE/CFSR arrays in `USflux_grids.inc` (`/canopy1/`, `/canopy2/`), which are also never populated — `ztan(30,1456,625)`, `zzan(30,1456,625)` and `htht(41,1456,625)` are ≈150 MB each of static arrays allocated for nothing. This is why the original Makefile build needed `-mcmodel=medium` / `-heap-arrays`. Removing these files from the link is the single biggest memory/startup win available in this codebase.

### 14.8 `getnetrad_simple` contains a self-overwriting statement
```fortran
rnetl = esfc*(esky-1.0)*sigma*(tak**4)
rnetl = esfc*esky*sigma*(tak**4) - esfc*sigma*((tak+4)**4)   ! overwrites the line above
```
The first assignment is dead. The second assumes surface temperature is `Ta + 4 K`.

### 14.9 `info/recipe/README.md` is a copy-paste stub
Its entire content is `# pyRTTOV` / "This project contains the python implementation of RTTOV." — leftover from a different conda recipe template. Ignore it.

### 14.10 `setup.py` builds `ALEXI_0.1`, the tree ships `ALEXI_0.2`
```python
version = "0.1"
mkPath  = os.path.join(process_dir, 'ALEXI_{}'.format(version))
```
Combined with `SConstruct.py`'s `include_path = .../ALEXI_0.2/src`, the packaging is internally inconsistent. The shipped binary matches `ALEXI_0.2` (it contains `landcover` and the free-form `USflux.f90`), so the 0.1 path in `setup.py` is stale.

### 14.11 Licence metadata contradicts itself
`setup.py`: `license='BSD 3-Clause'` but its classifiers include `License :: OSI Approved :: MIT License`. `meta.yaml` and `info/index.json`: `license: CCO`, `license_family: BSD`. `info/recipe/LICENSE.txt` is the authoritative text — check it before redistribution. Also a stale comment in `setup.py`: *"Uses dictionary comprehensions ==> 2.7 only"* under a Python 3.7 classifier.

### 14.12 Solver state carried across pixels
`hn0`, `psi0`, `fc0` (`common/initial/`) are `save`d and updated on every converged pixel, so each pixel's starting guess is the previous pixel's solution. Output ordering affects iteration counts and, at the margin, results. `ALEXI.f`'s `fc`-retry loop (`go to 1000`, `fcnew = fcnew + sign(0.01, 0.5-fcnew)`) is commented out, so a non-converged pixel is simply counted as `nfail`.

### 14.13 Hardcoded physics overrides worth knowing

| Location | Override |
|---|---|
| `USflux_run.f:97` | `iclass = 5` — land cover class from input is ignored |
| `USflux_run.f` | `iswater_inland = .FALSE.` — `ALEXI_water` unreachable |
| `USflux_run.f` | `fg = 1.0` — greenness fraction never varies |
| `USflux_run.f` | `emsoil = 0.94` fixed |
| `USflux_run.f` | 14 PBL heights hardcoded; `setPBLheights` (41 levels) unused |
| `ALEXI_utl.f:runinit` | `dtloc = 1.50` — the fc-dependent sunrise offset is computed then overwritten |
| `ALEXI_utl.f:set_constants` | `hn0 = 80 W/m²` initial guess |
| `USflux_cover.f` | `perennial = .FALSE.` forced |
| `USflux_cover.f` | `clump = 1.0` everywhere except `iclass == 12` (unreachable, since iclass is 5) |
| `landcover.f:canopyarch_lai` | `dispdh = 2/3`, `z0dh = 1/8` fixed |
| `USflux_utl.f:fao_PM` | reference-surface `albedo = 0.23` |
| `ALEXI_atmos.f:getresistance` | `val1 = 0.0025`, `val2 = 0.012` fixed (NDVI variants commented out) |

### 14.14 Unused / stale build artifacts in the source tree

`ALEXI_0.2/src/` ships 15 `.o` files and `.sconsign.dblite` from the original 2020 Linux build. Delete them before rebuilding; a stale `.sconsign.dblite` can make SCons skip recompilation.

### 14.15 `ALEXI_0.1` vs `ALEXI_0.2` divergence is undocumented
Both snapshots are present with the same file names. There is no changelog. Diff them before porting a fix, and pick one as canonical.

---

## 15. Corrections vs `README.md`

Audit of the pre-existing `README.md` against the source. Severity: **HIGH** = would misdirect work; **MED** = wrong detail; **LOW** = incomplete.

| # | Severity | `README.md` says | Source says |
|---|---|---|---|
| 1 | HIGH | "`USflux.f90` — **Alternative flux computation.** FAO Penman–Monteith formulation for direct ET estimation (non-inverse approach). Key subroutines: `usflux()`, `hourly_flux()`" | `USflux.f90` contains `subroutine USflux(MDATE,t,npoints,part)` and is **the main per-pixel driver loop** that dispatches `extract_input` → `ALEXI`/`ALEXI_water` → `clear_day_proc`/`cloudy_day_proc` → `screen_output`. There is no `usflux()` subroutine and no PM-based alternative solver. |
| 2 | HIGH | "`USflux_main.f` — Entry point. **Loops over all grid cells**, calls extract_input → ALEXI/USflux → output routines, tracks convergence stats" | `program main` does four `getarg`s, prints a banner, then `call landcover(...)` and `call USflux(...)`. The loop and the `nconv`/`nfail` counters live in `USflux.f90`. |
| 3 | HIGH | Execution-flow step 8: "`store_output_bin()` — Write results" | The only `call store_output_bin` is **commented out** (`USflux.f90:136`). Results are written by `screen_output` to a plain-text file. |
| 4 | HIGH | "`ALEXI_parm.inc` — Core model constants: von Kármán κ=0.4, π, cp=1010, PBL rise time, grid dimensions (ilg=1456, jlg=625), land cover classes" | `ALEXI_parm.inc` contains only `mli=41`, `ml=8000`, `mt=8000`, `BAD=-9999.0`, `SUCCESS=0`. κ/π/cp/`thrise` are assigned at runtime in `ALEXI_utl.f:set_constants`; grid dims and class counts are in `USflux_parm.inc`. |
| 5 | HIGH | "`USflux_grids.inc` — domain setup (CLAT=37.30°, CLON=−95.90°, dlat=0.04°)" | Those numbers appear **nowhere** in this tree. `CLAT, CLON, DX, DY, MINLAT…` are declared in `/canopy3/` and never assigned. |
| 6 | HIGH | Quick start: `conda install -c ashish-codebase projectmas_alexi` | The recipe channels are `conda-forge`, `bucricket`, `defaults`. `ashish-codebase` appears nowhere in the package metadata. |
| 7 | HIGH | Quick start: `pyalexi <date_YYYYDDD> <grid_name> <npoints> <part_number>` | `pyalexi` takes `--year/--doy/--grid_name` flags, not positional args — and it crashes before reaching the parser (§14.1). The positional form is `alexi_proc`'s interface. |
| 8 | HIGH | "Click any file link below to view its source code directly on GitHub" | There is no git repository here; `info/git` is a zero-byte file and `meta.yaml` points at a local path `/PycharmProjects/projectMASalexi`. |
| 9 | HIGH | Land-cover table of 14 classes attributed to `landcover.txt` | `landcover.txt` defines **7 usable classes** (Water, Forest, Shrubland, Grassland, Cropland, Bare_Ground, Urban) plus one malformed duplicate row. `nclass = 8`. The 14 MODIS-style names are hardcoded `tabdesc` strings in `USflux_cover.f` — and writing them overruns the array (§14.2). |
| 10 | MED | "`ALEXI.f` key subroutines: `alexi()`, `findflux_water()`, `growPBL()`" | `findflux_water()` lives in `ALEXI_water.f`. |
| 11 | MED | "`ALEXI_rad.f` key subroutines: `getradcomps()`, …" | `getradcomps()` lives in `USflux_rad.f`. |
| 12 | MED | FAO PM: `γ×900/(T+273)×u×(es−ea)` over `Δ + γ(1+0.34u)` | Code uses `γ·37·u·D/(Ta+273.15)` over `s + γ(1 + 0.24·u)` (hourly form), with `G = 0.1·Rn`. |
| 13 | MED | Relaxation factor "decreases from 0.25 to ~0.008" | `f = 0.25/(1+int(iter/20))`; the last value used before the 300-iteration bail is `0.25/15 ≈ 0.0167`. |
| 14 | MED | "`Ra = (log((z−d)/z0h) − ψm)/(κ·u*)`" presented as *the* Ra formula | That form is used only on the **bare-soil/water** branch (and uses `ψh`). The vegetated branch is `Ra = (xlog1−ψm)(xlog1−ψh)/(0.16·wind)`. |
| 15 | MED | `pbl_read.f` / `sfc_read.f` described as active readers | Both are compiled and linked but **never called** by the active pipeline. |
| 16 | MED | "`build.sh` — Shell script: runs SCons to compile Fortran" | `build.sh` runs `$PYTHON setup.py install`; `setup.py` is what invokes `scons`. |
| 17 | MED | "`info/recipe/SConstruct.py` — SCons build definition" | Actual path: `info/recipe/ALEXI_0.2/src/SConstruct.py`. |
| 18 | MED | Error-code table lists 0, 1, 2, 3, 9 | Code 9's message is commented out; code **30** (`Hn_new ≤ 0`) is missing from the table; 11/12/13 exist but are commented out at their set sites. |
| 19 | MED | "`landcover.f` key subroutines: `weighted_avg()`" | The entry point is `landcover(MDATE,tile,npoints,part)`; also `canopyarch_lai`, `loadclasstable`, `getdate`. |
| 20 | MED | Project structure shows `store/landcover.txt` at the top level and `info/meta.yaml` | Real paths: `info/recipe/ALEXI_0.2/store/landcover.txt` and `info/recipe/meta.yaml`. |
| 21 | LOW | "Dependencies: Python 3.7+" | Pinned `python >=3.7,<3.8.0a0`; `numpy 1.17.*`. |
| 22 | LOW | "~50 km resolution" | Not stated or derivable from the source; only `ilg=1456, jlg=625` is defined. |
| 23 | LOW | Omits `interp_LAI.f`, `SConstruct_testing.py`, `construct.yaml`, `conda_build_config.yaml`, `ALEXI_0.1/` | All present and relevant to a rebuild. |
| 24 | LOW | Omits that `info/recipe/README.md` is a `# pyRTTOV` stub | Confusing to any new contributor (§14.9). |
| 25 | LOW | Omits that `bin/alexi_proc` is a prebuilt **linux-64** ELF | Blocks any Windows/macOS run attempt. |

**Claims in `README.md` that check out:** the two-source energy-balance description; the inverse-for-`H` framing; `delhmx = 0.1 W/m²` and the 300-iteration cap; `Rs` and `Rx` formulas; `ALEXI_water` fixed albedo 0.1 / emissivity 0.99 / `G = 0.55·Rn`; `USflux_clear.f` / `USflux_cloud.f` / `USflux_cover.f` / `USflux_utl.f` / `ALEXI_utl.f` role summaries; `interp_LAI` linear-interpolation formula; `getradcomps` Cupid `RADIN4` lineage; `ALEXI_0.1` being the older snapshot.

---

## 16. Suggested next steps

Ordered by payoff per unit of effort.

1. **Decide the canonical source tree.** Delete or clearly quarantine `ALEXI_0.1/`, or document why both exist. Diff them once and record the delta.
2. **Delete the stale build artifacts** (`ALEXI_0.2/src/*.o`, `.sconsign.dblite`) and add a `.gitignore` for them.
3. **Fix `pyalexi`** so the CLI actually works: point the entry point at `pyalexi.pyalexi:main`, give `npoints` a real default, and pass `grid_name[0]` instead of the list.
4. **Fix `store/landcover.txt`**: renumber the trailing "Wetland" row to class 8 (or delete it), and decide whether class 5 is Cropland or Wetland.
5. **Fix the `tabdesc` overrun** in `USflux_cover.f` — either dimension `tabdesc` to 14 or drop assignments 9–14.
6. **Restore `iclass` from `pvars(5)`** in `extract_input` if you want per-pixel land cover to matter, and re-enable `iswater_inland` so `ALEXI_water` is reachable.
7. **Drop `pbl_read.f` and `sfc_read.f`** from the SCons source list (or drop the `/canopy1/`–`/canopy2/` arrays) to reclaim ~450 MB of static memory and remove the `-mcmodel` requirement.
8. **Re-enable a real output path.** Either uncomment `store_output_bin` (and create the binary grid plumbing it needs) or extend `screen_output`'s unit-103 line to include the daily totals, `xlec2`, `fpet`, and `iter`. Right now the only per-pixel product is five instantaneous values at t2.
9. **Add a smoke test**: a hand-built 3-pixel input set with known `Rn`, `SDN`, `TA`, `TRAD` and a two-level sounding, asserting convergence and a plausible `ET/Rn` ratio. Nothing in this tree is currently testable end-to-end.
10. **Rebuild on Linux/WSL** per §4 and confirm the rebuilt binary reproduces `bin/alexi_proc`'s behaviour before trusting any source edit.
11. **Reconcile licence metadata** (§14.11) against `LICENSE.txt` before publishing anywhere.
12. **Replace `info/recipe/README.md`** (the `# pyRTTOV` stub) with a pointer to this document.

---

## Appendix A — Quick fact sheet

| Item | Value | Source |
|---|---|---|
| Executable | `alexi_proc` (linux-64 ELF, 341,712 B) | `bin/`, `info/paths.json` |
| CLI signature | `alexi_proc YYYYDDD tile npoints part` | `USflux_main.f:86–95` |
| Grid | `ilg=1456`, `jlg=625` | `USflux_parm.inc` |
| Land-cover classes | `nclass=8` (7 populated) | `USflux_parm.inc`, `store/landcover.txt` |
| PBL input levels | `mli=41` declared; **14 used** | `ALEXI_parm.inc`, `USflux_run.f` |
| PBL interpolated levels | `ml=8000`; lookup `mt=8000` | `ALEXI_parm.inc` |
| Flux rise time | `thrise = 10800 s` | `ALEXI_utl.f:set_constants` |
| von Kármán κ | 0.4 | `ALEXI_utl.f:set_constants` |
| cp | 1010 J/kg·K | `ALEXI_utl.f:set_constants` |
| Bad-value sentinel | `BAD = -9999.0` | `ALEXI_parm.inc` |
| Convergence tolerance | `delhmx = 0.1 W/m²` | `ALEXI.f:22` |
| Max iterations | 300 | `ALEXI.f:64` |
| Initial guess | `hn0 = 80 W/m²`, `psi0 = 0` | `ALEXI_utl.f:set_constants` |
| Blending height `zta` | 50 m | `USflux_run.f` |
| Wind measurement height `refhtw` | 30 m → reset to 50 m | `USflux_run.f`, `ALEXI_utl.f:runinit` |
| Wind clamp | 3–20 m/s | `ALEXI_utl.f:runinit` |
| Sunrise offset `dtloc` | 1.50 h (hardcoded) | `ALEXI_utl.f:runinit` |
| Leaf / soil emissivity | 0.97 / 0.94 | `ALEXI_utl.f:runinit`, `USflux_run.f` |
| Hourly loop length | `nohr = 24` | `USflux_parm.inc` |
| Daily flux units | MJ/m²/day | `USflux_utl.f:integrate` |
| Instantaneous flux units | W/m² | throughout |

## Edit History

- **2025 — v2 created.** Full source audit of `info/recipe/ALEXI_0.2/src/` (15 Fortran units, 9 include files), the conda packaging layer, and the `pyalexi` wrapper. Written as the handoff baseline; §15 records every divergence found in the original `README.md`.
