GRMONTY: A relativistic Monte Carlo code
==========================================

Based on [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D). Originally downloaded from [Astrophysical Code Library](http://rainman.astro.illinois.edu/codelib/) @ [UI](http://illinois.edu).

GRMONTY is parallelized using [OpenMP](https://en.wikipedia.org/wiki/OpenMP).

This version is configured to use input from [`harm2d`](http://rainman.astro.illinois.edu/codelib/codes/ham2d/src/), which is also available on this site.

# quick start

unpack the tarball:

    tar -xzvf grmonty.tgz

make (requires openmp enabled gcc):

    make

set number of threads (example is for csh and 8 threads):

    setenv OMP_NUM_THREADS 8

run the code on the supplied harm output file:

    ./grmonty 5000000 dump019 4.e19 >&! grmonty.err

Arguments are:

- estimate of photon number (actual number is probabilistic due to scattering)
- harm dump file for model
- mass unit (few x 10^19 is appropriate for Sgr A*)

This will output spectrum to `grmonty.spec`  which should be identical to `grmonty_spec_verify`.

Then use SM script `plspec` to plot up broad-band spectral energy distribution.

To calculate spectra from other sources, replace `harm_model.c`
with your own source model.  Begin by modifying `harm_model.c`.

You must supply

```
init_model 
make_super_photon
bias_func
get_fluid_params
report_spectrum
stop_criterion
record_criterion

gcon_func 
gcov_func 
get_connection
```

in the model file.

# TODO

- [ ] make it work with HARMPI/HARM3D
- [ ] make it work with RAISHIN
- [ ] add LICENSE

# References

- [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D)
- TBD

# LICENSE 

`grmonty` is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the `LICENSE` file or the [GNU General Public License](http://www.gnu.org/licenses/) for more details.
