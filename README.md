# ALEXI (Atmospheric Land Exchange Inverse) — Source Code Reference

A Conda-packaged implementation of the **ALEXI model** for partitioning satellite-observed land surface energy fluxes into evapotranspiration (latent heat), sensible heat, and soil heat flux. Developed by [Martha Anderson](https://ams.org/ams-blog/2017/03/14/martha-anderson) at the US Agricultural Research Service (USDA-ARS).

> **Click any file link below** to view its source code directly on GitHub. Use this as a navigational map through the ~25 Fortran source files that make up the model.

---

## Table of Contents

- [What ALEXI Does](#what-alexi-does)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Architecture & Data Flow](#architecture--data-flow)
- [Execution Flow](#execution-flow)
- [File Reference — Include Files](#file-reference---include-files)
- [File Reference — Main Execution Flow](#file-reference---main-execution-flow)
- [File Reference — Core Physics](#file-reference---core-physics)
- [File Reference — Utilities & Data I/O](#file-reference---utilities--data-io)
- [Key Equations](#key-equations)
- [Build System](#build-system)

---

## What ALEXI Does

ALEXI is a **two-source energy balance model** that estimates real-time evapotranspiration (ET) over large areas using satellite observations. It works by:

1. Taking **GOES radiometric surface temperatures** at two times (~10:30 and ~13:30 local time)
2. Using **meteorological profiles** (air temp, humidity, wind, pressure) from reanalysis data
3. Solving an **inverse problem**: finding the sensible heat flux (H) that is consistent with the observed temperature difference between the two times and the growth of the planetary boundary layer (PBL)
4. Partitioning net radiation into **latent heat (LE)**, **sensible heat (H)**, and **soil heat flux (G)**

The model runs on a continental US domain at ~50 km resolution (1456 × 625 grid cells).

---

## Quick Start

```bash
# Install from Conda
conda install -c ashish-codebase projectmas_alexi

# Run the Python CLI wrapper
pyalexi <date_YYYYDDD> <grid_name> <npoints> <part_number>
```

The Python entry point [`bin/pyalexi`](bin/pyalexi) calls the compiled Fortran binary [`bin/alexi_proc`](bin/alexi_proc) via `subprocess`.

---

## Project Structure

```
projectmas_alexi-0.1.0/
├── README.md                          ← You are here
├── bin/
│   ├── pyalexi                        ← Python CLI wrapper
│   └── alexi_proc                     ← Compiled Fortran binary (ELF)
├── info/
│   ├── meta.yaml                      ← Conda package metadata
│   ├── index.json                     ← Package index
│   └── recipe/
│       ├── ALEXI_0.1/src/             ← Older version source
│       ├── ALEXI_0.2/src/             ← Active version source (see below)
│       ├── pyalexi/                   ← Python package source
│       │   ├── __init__.py
│       │   └── pyalexi.py             ← CLI entry point
│       ├── build.sh                   ← Conda build script
│       ├── SConstruct*                ← SCons build files
│       └── meta.yaml.template         ← Conda template
├── lib/                               ← Installed Python packages
└── store/                             ← Land cover lookup tables
    └── landcover.txt                  ← 14 land cover class parameters
```

**All Fortran source lives in:** [`info/recipe/ALEXI_0.2/src/`](info/recipe/ALEXI_0.2/src/)

---

## Architecture & Data Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│                        INPUT DATA SOURCES                            │
│  CFSR/WRF: TA, EA, PRES, WIND, SDN, XLWDN (hourly)                 │
│  GOES:     Radiometric temperatures (TRAD1, TRAD2) at ~10:30/13:30  │
│  Satellite: LAI, fraction cover, land cover class                   │
│  Sounding: PBL temperature/pressure profiles                        │
└───────────────────────────────┬─────────────────────────────────────┘
                                │
                    ┌───────────▼───────────┐
                    │   USflux_main.f       │
                    │   (Driver Program)    │
                    │   Loop over grid cells│
                    └───────────┬───────────┘
                                │
              ┌─────────────────┼─────────────────┐
              │                 │                 │
      ┌───────▼───────┐ ┌──────▼───────┐ ┌───────▼───────┐
      │ USflux_run.f  │ │ ALEXI.f      │ │ USflux.f90    │
      │ extract_input │ │ Inverse model│ │ PM-based      │
      │ getdgmt       │ │ iterative Hn │ │ direct solve  │
      │ screen_output │ │ solver       │ │               │
      │ store_output  │ └───────┬───────┘ └───────┬───────┘
      └───────────────┘       │                 │
              │          ┌─────▼─────┐    ┌─────▼─────┐
              │          │ ALEXI_atmos│   │ USflux_utl │
              │          │ .f         │    │ .f         │
              │          │ (Aerodynamic│    │ (Hourly    │
              │          │  Resistance)│    │  fluxes)   │
              │          └─────┬───────┘    └─────┬───────┘
              │                │                  │
              │         ┌──────▼──────┐   ┌──────▼──────┐
              │         │ ALEXI_rad.f │   │ USflux_clear│
              │         │ (Radiation) │   │ .f / _cloud  │
              │         └─────────────┘   │ .f           │
              │                           └──────┬───────┘
              │          ┌───────────────────────┤
              ▼          ▼                       ▼
      ┌─────────────────────────────────────────────────┐
      │              OUTPUT FILES                        │
      │  Binary flux matrices: Rn, LE, H, G (daily)     │
      │  Text diagnostics per grid cell                 │
      └─────────────────────────────────────────────────┘
```

---

## Execution Flow

Each grid cell goes through this pipeline:

1. **[`USflux_main.f`](info/recipe/ALEXI_0.2/src/USflux_main.f)** — Driver loops over all (ia, ja) grid cells
2. **[`USflux_run.f`](info/recipe/ALEXI_0.2/src/USflux_run.f)** → `extract_input()` — Reads raw input for one pixel
3. **[`ALEXI.f`](info/recipe/ALEXI_0.2/src/ALEXI.f)** or **[`USflux.f90`](info/recipe/ALEXI_0.2/src/USflux.f90)** — Core physics solver
4. **[`ALEXI_atmos.f`](info/recipe/ALEXI_0.2/src/ALEXI_atmos.f)** — Aerodynamic resistance calculations
5. **[`ALEXI_rad.f`](info/recipe/ALEXI_0.2/src/ALEXI_rad.f)** — Radiation partitioning
6. **[`USflux_utl.f`](info/recipe/ALEXI_0.2/src/USflux_utl.f)** → `hourly_flux()` — Hourly energy balance
7. **[`USflux_clear.f`](info/recipe/ALEXI_0.2/src/USflux_clear.f)** / **[`USflux_cloud.f`](info/recipe/ALEXI_0.2/src/USflux_cloud.f)** — Daily flux aggregation
8. **[`USflux_run.f`](info/recipe/ALEXI_0.2/src/USflux_run.f)** → `store_output_bin()` — Write results

---

## File Reference — Include Files

Include files define shared data structures (common blocks, parameters) used across multiple source files. Click any link to view:

| File | Purpose |
|------|---------|
| [`ALEXI_parm.inc`](info/recipe/ALEXI_0.2/src/ALEXI_parm.inc) | Core model constants: von Karman constant (κ=0.4), π, cp=1010 J/kg·K, PBL rise time (3 hrs), grid dimensions (ilg=1456, jlg=625), land cover classes |
| [`ALEXI_f90.inc`](info/recipe/ALEXI_0.2/src/ALEXI_f90.inc) | Fortran 90 type/parameter declarations |
| [`ALEXI.inc`](info/recipe/ALEXI_0.2/src/ALEXI.inc) | Common blocks for core ALEXI variables: canopy, data1/data2 (t1/t2 observation times), PBL, flags |
| [`USflux_parm.inc`](info/recipe/ALEXI_0.2/src/USflux_parm.inc) | US Flux parameters: grid sizes, land cover class definitions, NCAR graphics params |
| [`USflux_grids.inc`](info/recipe/ALEXI_0.2/src/USflux_grids.inc) | 2D grid arrays for all model fields; domain setup (CLAT=37.30°, CLON=-95.90°, dlat=0.04°) |
| [`USflux_dir.inc`](info/recipe/ALEXI_0.2/src/USflux_dir.inc) | Directory paths: HOMEDIR, INDIR (inputs), TABDIR (landcover tables), OUTDIR (fluxes) |
| [`USflux.inc`](info/recipe/ALEXI_0.2/src/USflux.inc) | US Flux common blocks: daily fluxes, hourly data, observation fluxes, moisture stress flags |
| [`date.inc`](info/recipe/ALEXI_0.2/src/date.inc) | Date variables in common `/canopy0/`: JDATE, MDATE (YYYYDDD format), CDATE, YDATE |
| [`flag.inc`](info/recipe/ALEXI_0.2/src/flag.inc) | Flag array for gridded processing flags |

---

## File Reference — Main Execution Flow

These files control the overall program execution and data handling.

### Driver & Orchestration

| File | Role |
|------|------|
| [`USflux_main.f`](info/recipe/ALEXI_0.2/src/USflux_main.f) | **Entry point.** Loops over all grid cells, calls extract_input → ALEXI/USflux → output routines, tracks convergence stats (nconv, nfail) |
| [`USflux_run.f`](info/recipe/ALEXI_0.2/src/USflux_run.f) | **Input extraction & output.** `extract_input()` reads raw data per pixel; `screen_output()` writes diagnostics; `store_output_bin()` writes binary results |

### Python Entry Point

| File | Role |
|------|------|
| [`pyalexi/pyalexi.py`](info/recipe/ALEXI_0.2/src/pyalexi/pyalexi.py) | Python CLI wrapper that uses `subprocess` to call the compiled Fortran binary `alexi_proc` |
| [`bin/pyalexi`](bin/pyalexi) | Installed Python CLI entry point script |

---

## File Reference — Core Physics

These files contain the mathematical heart of the model.

### ALEXI Inverse Model (Fortran 77)

| File | Role | Key Subroutines |
|------|------|-----------------|
| [`ALEXI.f`](info/recipe/ALEXI_0.2/src/ALEXI.f) | **Core inverse solver.** Iteratively finds sensible heat flux (Hn) using adaptive relaxation until convergence (ΔH < 0.1 W/m²). Max 300 iterations. | `alexi()`, `findflux_water()`, `growPBL()` |
| [`ALEXI_atmos.f`](info/recipe/ALEXI_0.2/src/ALEXI_atmos.f) | **Aerodynamic resistance & PBL.** Computes Ra, Rs, Rx using Monin-Obukhov stability theory (Kader & Yaglom 1990). Canopy architecture from Massman & Weil (1996). | `getresistance()`, `psimhn_ky()`, `canopyarch()` |
| [`ALEXI_rad.f`](info/recipe/ALEXI_0.2/src/ALEXI_rad.f) | **Radiation partitioning.** Shortwave (visible/NIR, direct/diffuse) and longwave fluxes using Goudriaan 1988 canopy reflectance/transmission models. Brutsaert equation for downward longwave. | `getradcomps()`, `getradprops()`, `getnetrad()`, `getnetradnight()`, `getxlwdn()` |
| [`ALEXI_water.f`](info/recipe/ALEXI_0.2/src/ALEXI_water.f) | **Water body solver.** Specialized physics for lakes/oceans: fixed albedo=0.1, emissivity=0.99, soil heat = 55% of net radiation. | `alexi_water()`, `findflux_water()`, `getnetrad_water()` |

### US Flux Model (Penman-Monteith Based)

| File | Role | Key Subroutines |
|------|------|-----------------|
| [`USflux.f90`](info/recipe/ALEXI_0.2/src/USflux.f90) | **Alternative flux computation.** FAO Penman-Monteith formulation for direct ET estimation (non-inverse approach). | `usflux()`, `hourly_flux()` |
| [`USflux_clear.f`](info/recipe/ALEXI_0.2/src/USflux_clear.f) | **Clear-sky daily flux aggregation.** Computes system, soil, and canopy fluxes from hourly totals. FAO Penman-Monteith reference ET for stress factors. | `clear_day_proc()`, `daily_flux_clear_SDN()` |
| [`USflux_cloud.f`](info/recipe/ALEXI_0.2/src/USflux_cloud.f) | **Cloudy day handling.** Sets all daily outputs to BAD (gap-filled in post-processing). | `cloudy_day_proc()`, `daily_flux_cloudy()` |
| [`USflux_cover.f`](info/recipe/ALEXI_0.2/src/USflux_cover.f) | **Vegetation properties.** Land cover classification, LAI-to-cover conversion, canopy architecture parameters, clumping factors. Reads `landcover.txt`. | `load_tables()`, `cover_props()`, `getfveg()` |
| [`USflux_rad.f`](info/recipe/ALEXI_0.2/src/USflux_rad.f) | **US Flux radiation utilities.** Cupid RADIN4 strategy for SDN partitioning; simplified net radiation from solar + air temp. | `getradcomps()`, `getnetrad_simple()` |
| [`USflux_utl.f`](info/recipe/ALEXI_0.2/src/USflux_utl.f) | **Hourly flux computation & integration.** 24-hour loop: validates inputs, computes solar zenith, partitions radiation, estimates soil heat, integrates to daily totals. FAO Penman-Monteith reference ET. | `hourly_flux()`, `fao_PM()`, `integrate()`, `getTdepart2()` |
| [`landcover.f`](info/recipe/ALEXI_0.2/src/landcover.f) | **Land cover classification.** Reads vegetation input files, calculates LAI, fractional cover, calls weighted_avg. | `weighted_avg()` |
| [`interp_LAI.f`](info/recipe/ALEXI_0.2/src/interp_LAI.f) | **LAI interpolation program.** Standalone: linearly interpolates between two LAI datasets (val = v1 + f × (v2 - v1)). | `INTERP_LAI` (main program) |

---

## File Reference — Utilities & Data I/O

| File | Role | Key Subroutines |
|------|------|-----------------|
| [`ALEXI_utl.f`](info/recipe/ALEXI_0.2/src/ALEXI_utl.f) | **Utilities & initialization.** Constants (κ, π, cp, thrise), PBL profile interpolation, input validation with range checks, error code mapping. | `set_constants()`, `getprofile()`, `runinit()`, `checkinput()`, `checkvalue()` |
| [`pbl_read.f`](info/recipe/ALEXI_0.2/src/pbl_read.f) | **PBL profile reader.** Reads CFSR temperature/pressure profiles, converts potential temp to actual via hypsometric equation, interpolates to 200m levels. | `pbl_read()`, `fill_hgt_domain()` |
| [`sfc_read.f`](info/recipe/ALEXI_0.2/src/sfc_read.f) | **Surface data reader.** Reads gridded CFSR/WRF surface meteorology (TA, humidity, pressure, wind, longwave). Sets up domain lat/lon arrays. | `sfc_read()`, `arraynav()` |

---

## Key Equations

### Energy Balance
```
Rn = LE + H + G
```
Net radiation = latent heat + sensible heat + soil heat flux

### Inverse Model (ALEXI.f)
```
Adaptive relaxation:  hn(new) = hn(old) + f × (hn(new) - hn(old))
where f = 0.25 / (1 + int(iter/20))   [decreases from 0.25 to ~0.008]
Convergence: |ΔH| < 0.1 W/m² or iter > 300
```

### Aerodynamic Resistance (ALEXI_atmos.f)
```
Ra = (log((z-d)/z0h) - ψm) / (κ × u*)
Rs = 1 / (a1 × ΔT^0.33 + a2 × us)
Rx = 180 × sqrt(xl/udz) / LAI
```

### Radiation (ALEXI_rad.f — Goudriaan 1988)
```
Canopy albedo = f(visible, NIR) × [direct + diffuse] reflectance
Transmittance = exp(-k × clump × LAI)
Net radiation = ε_sky × σ × T_sky⁴ - ε_surf × σ × T_surf⁴ + SW_in × (1 - α)
```

### FAO Penman-Monteith (USflux_utl.f)
```
λE = [0.408×Δ×(Rn-G) + γ×900/(T+273)×u×(es-ea)] / [Δ + γ×(1+0.34×u)]
```

---

## Build System

The project uses **conda-build** with **SCons** for Fortran compilation:

| File | Purpose |
|------|---------|
| [`info/recipe/meta.yaml`](info/recipe/meta.yaml) | Conda package metadata (name, version, dependencies, build script) |
| [`info/recipe/build.sh`](info/recipe/build.sh) | Shell script: runs SCons to compile Fortran → `alexi_proc` |
| [`info/recipe/SConstruct.py`](info/recipe/ALEXI_0.2/src/SConstruct.py) | SCons build definition: compiles all `.f` files with gfortran |

**Dependencies:** Python 3.7+, NumPy 1.17, GCC/GFortran 9.3+

---

## Land Cover Classes (from `landcover.txt`)

| Class | Description |
|-------|-------------|
| 1 | Water (and Goodes Interrupted Space) |
| 2 | Evergreen Needleleaf Forest |
| 3 | Evergreen Broadleaf Forest |
| 4 | Deciduous Needleleaf Forest |
| 5 | Deciduous Broadleaf Forest |
| 6 | Mixed Cover |
| 7 | Woodland |
| 8 | Wooded Grassland |
| 9 | Closed Shrubland |
| 10 | Open Shrubland |
| 11 | Grassland |
| 12 | Cropland |
| 13 | Bare Ground |
| 14 | Urban and Built-Up |

---

## Error Codes (ALEXI.f)

| Code | Meaning |
|------|---------|
| 0 | **CONVERGED** — successful solution |
| 1 | Hn did not converge within 300 iterations |
| 2 | Bad input data — bail |
| 3 | Exceeded max PBL profile layer |
| 9 | Net radiation < 0 — physically impossible |

---

*Generated from analysis of ALEXI v0.1.0 source code. For questions or contributions, open an issue on this repository.*
