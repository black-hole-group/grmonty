# AGENTS.md - GRMONTY Project Guidelines

## Project Overview

GRMONTY is a relativistic Monte Carlo radiative transfer code for computing spectra from GRMHD simulations near black holes. It calculates emergent spectra via thermal synchrotron and inverse Compton scattering.

Reference: Dolence et al. 2009, ApJS, 184, 387

## Directory Structure

```
grmonty/
├── src/
│   ├── main.c              # Entry point
│   ├── physics/            # Radiation physics
│   ├── geometry/           # Geodesics, tetrads
│   └── model/              # HARM model implementation
├── include/                # Header files
├── data/                   # Input data files
├── test/                   # Test reference files
├── scripts/                # Utility scripts
├── build/                  # Build artifacts (generated)
└── makefile
```

## Build/Lint/Test Commands

### Build
```bash
make                  # Build the project (requires OpenMP-enabled gcc)
make clean            # Remove build artifacts
```

### Running
```bash
export OMP_NUM_THREADS=8           # Set threads (bash)
./grmonty 5000000 data/dump019 4.e19    # Run: photon_estimate, harm_dump_file, mass_unit
```

### Verification
```bash
diff grmonty.spec test/grmonty_spec_verify   # Compare output to reference
```

### Linting
```bash
make clean && make 2>&1 | head -50       # Check for warnings (-Wall enabled)
```

## Code Style Guidelines

### File Structure
- All `.c` files include `<decs.h>` as first header (angle brackets for include/)
- Model-specific files also include `<harm_model.h>`
- Headers: `decs.h` (main), `constants.h` (CGS physical constants), `harm_model.h` (model-specific)

### Imports
```c
#include <decs.h>        // Always first
#include <harm_model.h>  // For model-specific files (after decs.h)
```
Standard library includes (`<stdio.h>`, `<stdlib.h>`, `<math.h>`, GSL headers, `<omp.h>`) are in `decs.h`.

### Formatting
- Indent with tabs
- K&R style braces, opening brace on same line
- Maximum line length: ~80 characters

### Naming Conventions
| Type | Convention | Example |
|------|------------|---------|
| Functions | lowercase_underscores | `track_super_photon()`, `get_fluid_nu()` |
| Variables | lowercase_underscores | `thetae`, `nu`, `theta` |
| Constants/Macros | UPPERCASE_UNDERSCORES | `NDIM`, `SMALL`, `N_EBINS` |
| Structures | `struct of_<name>` | `struct of_photon`, `struct of_geom` |
| Globals | lowercase, extern in decs.h | `extern double a;` |
| Thread-local | Use `#pragma omp threadprivate` | `gsl_rng *r; #pragma omp threadprivate(r)` |

### Types
- `double` for all floating-point physics calculations
- `int` for indices and counters
- Loop indices: `i`, `j`, `k`, `l` convention

### Useful Macros
```c
#define DLOOP for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++)
#define INDEX(i,j,k) (NPRIM*((k) + N3*((j) + N2*(i))))
#define SMALL 1.e-40   // Prevent division by zero
```

### Error Handling
```c
if (fp == NULL) {
    fprintf(stderr, "can't open sim data file\n");
    exit(1);
}

if (isnan(ph->w) || isnan(ph->E)) {
    fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
    return;
}
```

### OpenMP Parallelization
```c
#pragma omp parallel private(ph)
{
    while (1) {
#pragma omp critical (MAKE_SPHOT)
        { make_super_photon(&ph, &quit_flag); }
        // ...
    }
}

#pragma omp atomic
N_superph_made += 1;
```

### Numerical Safety
- Use `SMALL` to prevent division by zero: `x / (y + SMALL)`
- Check bounds before array access
- Use `fabs()` for absolute values with doubles
- Use `isnan()` and `isinf()` checks

### Physical Constants
Defined in `include/constants.h` in CGS units: `EE` (electron charge), `CL` (speed of light), `ME` (electron mass), `MP` (proton mass), `HPL` (Planck constant), `KBOL` (Boltzmann constant).

### GSL Usage
- Random number generation: `gsl_rng_mt19937` (Mersenne Twister)
- Random direction: `gsl_ran_dir_3d(r, &x, &y, &z)`
- Use GSL for special math: Bessel functions, integration, chi-squared distribution

### Adding New Physics
1. Add function prototype to `include/decs.h`
2. Implement in new `.c` file in appropriate subdirectory
3. Add to `SRCS` and create build rule in makefile
4. Follow existing patterns for lookup tables

### Memory Management
- Use `INDEX` macro for 3D array access
- Use `malloc_rank1`, `malloc_rank2`, `malloc_rank2_cont` helpers in `src/model/harm_utils.c`
- No dynamic allocation in inner loops (performance critical)

### Output
- Spectrum output: `grmonty.spec`
- Progress to stderr with `fprintf(stderr, ...)`
- Use `fflush(stderr)` after progress messages

### Code Header
All source files must include the GPL copyright header (see existing files).
