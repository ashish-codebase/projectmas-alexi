# ALEXI C# port

C# (C# 2024 / .NET 8) port of the ALEXI_0.1 Fortran model
(`info/recipe/ALEXI_0.1/src/`). Pure .NET — no native dependencies — so it
opens and runs in **Visual Studio Community Edition**.

## Opening in Visual Studio

Open `Alexi.sln` (legacy format; `Alexi.slnx` is equivalent). Requirements:

- Visual Studio 2022 (17.x) or newer, or Visual Studio Code with the .NET 8 SDK.
- No native toolchain or gfortran needed.

Debug the CLI from `src/Alexi.Cli` with arguments (see below).

## Command line

```
alexi <date_YYYYDDD> <grid_tag> <npoints> <part>
      [--input-dir D] [--output-dir D] [--table-dir D]
      [--log-dir D] [--albedo-dir D] [--lai-dir D] [--store-binary]
```

- `date_YYYYDDD` — e.g. `2008172` (2008, day 172).
- `grid_tag`, `part` — appended to input filenames exactly as the Fortran driver does.
- `npoints` — number of pixel records to process.
- Defaults: `--input-dir ./INPUTS/ALEXI_INPUT/`, `--output-dir ./OUTPUTS/ALEXI_OUTPUT/`,
  `--table-dir info/recipe/ALEXI_0.1/`, `--log-dir ./ALEXI_LOG`.
- `--store-binary` enables `store_output_bin` (commented out in the shipped Fortran; off by default).

Subcommands (offline utilities):

```
alexi interp-lai <inputDir>
```

Example:

```
dotnet run --project src/Alexi.Cli -- 2008172 _npt 100 us
```

## Build & test

```
dotnet build Alexi.sln
dotnet test          # Alexi.Core.Tests + Alexi.Golden.Tests
```

## Layout

| Project | Contents |
| --- | --- |
| `src/Alexi.Core` | Model physics: atmos, rad, water, two-source solver, cover, clear/cloudy flux, utilities, state records |
| `src/Alexi.IO` | Driver (`USflux.f90`), input/screen/binary I/O (`USflux_run.f`), landcover stage, offline readers (`pbl_read`, `sfc_read`) |
| `src/Alexi.Cli` | `USflux_main.f` front end |
| `tests/*` | Hand-derived unit tests (NOT Fortran parity — no oracle on this host) |

## Caveats

- Tests are hand-derived expectations, not golden output from the Fortran binary
  (`bin/alexi_proc` is a Linux ELF and cannot run on this Windows host).
- Known Fortran defects are preserved verbatim and documented in file headers
  (see `handoff.md` for the full list).
