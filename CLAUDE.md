# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Run

```bash
make                    # build -> produces ./grmonty binary (requires OpenMP-enabled gcc)
make clean              # remove object files
make clean && make 2>&1 | head -50  # check for compiler warnings
```

Run with sample data:
```bash
export OMP_NUM_THREADS=8
./grmonty 5000000 data/dump019 4.e19
# Arguments: photon_estimate  harm_dump_file  mass_unit(CGS)
```

Verify output:
```bash
diff grmonty.spec test/grmonty_spec_verify
```

Plot spectrum with Python:
```python
import nmmn.sed
s = nmmn.sed.SED()
s.grmonty('grmonty.spec')
plot(s.lognu, s.ll)
```

## Architecture

GRMONTY is a relativistic Monte Carlo radiative transfer code (Dolence et al. 2009) that computes emergent spectra from GRMHD simulation dumps via thermal synchrotron emission and inverse Compton scattering.

**Data flow:**
1. `init_harm_data.c` reads a HARM binary dump → primitive MHD variables `p[var][i][j]`
2. `init_geometry.c` precomputes `gcov`, `gcon`, `sqrtg` at every grid zone
3. `harm_model.c` sets physical units and builds emissivity/weight lookup tables
4. OpenMP parallel region (one thread per core):
   - `make_super_photon()` samples emission location, frequency, and weight
   - `track_super_photon()` transports each photon along geodesics, computing optical depths and triggering Compton scatter events
   - `record_super_photon()` bins escaped photons into `spect[N_THBINS=6][N_EBINS=200]`
5. Spectra are reduced across threads and written to `grmonty.spec` (νLν in solar luminosities, binned by viewing angle and log energy)

**Key design patterns:**
- **Super-photons**: weighted Monte Carlo packets; weight decays as `exp(-τ_abs)`. Russian roulette kills packets below `WEIGHT_MIN` or boosts them by `ROULETTE=1e4`.
- **Scattering bias**: `bias_func` artificially increases scattering rate with weight compensation to improve inverse Compton sampling.
- **Pluggable model layer**: the transport core is model-independent. To add a new plasma source, replace `harm_model.c` supplying the same ~12 function signatures (`init_model`, `make_super_photon`, `bias_func`, `get_fluid_params`, `stop_criterion`, `record_criterion`, `report_spectrum`, `gcon_func`, `gcov_func`, `get_connection`, etc.).
- **MKS coordinates**: the metric and Christoffel symbols are implemented analytically for a Kerr BH in Modified Kerr-Schild coordinates used by HARM.
- **Compile-time configuration**: frequency range (`NUMIN`/`NUMAX`), temperature limits, `TP_OVER_TE`, grid resolution (`N_EBINS`, `N_THBINS`), and `RMAX` are `#define` constants in `decs.h`.

**Key files:**
| File | Role |
|------|------|
| `src/main.c` | `main()`: argument parsing, OMP parallel loop, spectrum reduction |
| `include/decs.h` | Master header: all structs, extern globals, macros, function prototypes |
| `include/constants.h` | CGS physical constants |
| `include/harm_model.h` | Model-specific declarations |
| `src/model/harm_model.c` | Pluggable model layer (HARM-specific metric + fluid interface) |
| `src/model/track_super_photon.c` | Core MC transport: geodesic stepping, opacity integration, scattering |
| `src/geometry/geodesics.c` | Predictor-corrector geodesic integrator (`push_photon`) |
| `src/geometry/tetrads.c` | Tetrad algebra for frame transformations |
| `src/physics/compton.c` | Klein-Nishina sampling, Mersenne Twister RNG via GSL |
| `src/physics/scatter_super_photon.c` | Compton event: boost to electron rest frame, scatter, boost back |
| `src/physics/radiation.c` | Model-independent opacity/emissivity utilities |
| `src/physics/jnu_mixed.c` | Thermal synchrotron emissivity + lookup tables (uses GSL Bessel functions) |
| `src/physics/hotcross.c` | Hot Compton cross-section 2D lookup table; reads `data/hotcross.dat` |
| `src/model/harm_utils.c` | Grid helpers: storage allocation, coordinate mapping, bilinear interpolation |

**Dependencies:** GSL (RNG, Bessel functions, integration), OpenMP (`-fopenmp`), libm. No MPI or GPU code in this repo (see [GPUmonty](https://github.com/black-hole-group/gpumonty) for GPU port).

## Code Style

- K&R C style; tabs for indentation; ~80 character line limit
- All `.c` files include `<decs.h>` first (angle brackets), then `<harm_model.h>` if model-specific
- Standard library headers live in `include/decs.h`; do not re-include them in `.c` files
- Functions and variables: `lowercase_with_underscores`
- Macros and constants: `UPPERCASE_WITH_UNDERSCORES`
- Structs: `struct of_<descriptive_name>`
- Globals defined in `grmonty.c`, declared `extern` in `decs.h`
- Thread-local globals: `#pragma omp threadprivate(var)`
- Use `SMALL` (1.e-40) to guard against division by zero
- Floating point: `double` throughout; `fabs()` for absolute values
- Fatal errors: `fprintf(stderr, ...); exit(1);`
- Progress/diagnostics: `fprintf(stderr, ...); fflush(stderr);`
- New source files go in the appropriate subdirectory (`src/physics/`, `src/geometry/`, or `src/model/`), must be added to `SRCS` in `makefile`, and prototypes added to `include/decs.h`
- All source files carry the GPL v3 copyright header
