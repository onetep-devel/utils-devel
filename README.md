# Utilities
Utility scripts and other resources for the ONETEP linear-scaling DFT code.

## Interpretation of ONETEP runs

- ### `summarise`: Makes a summary of an output produced with `output_detail VERBOSE`.

  Bash/awk script by Nicholas Hine. Extracts the results of the NGWF CG optimisation steps from an output file 
  (which may still be running) and outputs them in a format as if you were 
  running with `output_detail BRIEF` or looking at the calculation summary.

- ### `geomconv`: Helps investigate convergence indicators of ONETEP geometry optimisation.

  Bash/awk script by Nicholas Hine. Extracts the convergence indicators of a ONETEP BFGS 
  geometry optimisation calculation from the output file, and compares them to the convergence 
  tolerances. The results are coloured to indicate which parameters are converged and 
  which are not.

  *Note*: This script needs another scripts to be present in your `PATH`: `summarise`.

- ### `kern_vs_ngwf`: Helps investigate convergence indicators of a ONETEP single point energy run.

  Bash/awk script by Nicholas Hine. 
  Views the convergence of an ongoing ONETEP calculation, and displays various coloured
  indicators as to the convergence of the calculation. Examine the results to view
  decreases in energy, commutator, and ratio of optimisation energy gain from kernel and 
  NGWF optimisation and spot common problems (requires some experience to use successfully).


## Format conversion

- ### `dat2xsf`: Converts `.dat` files for ONETEP into `.xsf` files for direct loading into xcrysden.

  Bash/awk script by Nicholas Hine. See the header of this script for detailed instructions.

- ### `geom2xsf`: Converts `.geom` or `.md` files produced by ONETEP into `.xsf` files for loading into xcrysden.

  Bash/awk script by Nicholas Hine. See the header of this script for detailed instructions.

- ### `geom2xyz`: Converts `.geom` or `.md` files from ONETEP or CASTEP into `.xyz` for animation in xmol.

  Perl script by Greg Pearce. This is from 2002 and no longer maintained. It might not function correctly in this day and age.

- ### `usp2upf`: Converts a pseudopotential in the CASTEP `.usp` form to unified pseudopotential format.

  Fortran program by PWSCF group. Distributed under the terms of the GNU General Public license.

## Molecular dynamics

- ### `md2thermo`: Calculates thermodynamical quantities from a ONETEP .md file.

  Bash/awk script by Jacek Dziedzic. Can be used to process an MD trajectory (`.md` file) to get statistical information on energy conservation and temperature.

  *Note*: This script needs several other scripts to be present in your `PATH`. These are: `getcol`, `getline`, `findmin`, `findmax`, `findavg`, `findrmsavg`, `findabsmax`, `findstddev`, `format_as`, and `linear_fit`. These scripts are also provided here.

- ### `md_clean_duplicates`: Cleans up an `.md` trajectory from multiple restarts.

  Bash/awk script by Jacek Dziedzic. Reads an `.md` trajectory file, copying input to output, except for 
  duplicate steps (those with identical times), which are removed. Only the first instance of any time is retained. 
  This is useful for processing `.md` files combined from multiple runs, where global history has been used, 
  leading to short time intervals present at the end of one run, and then again at the beginning of a restared run.

## Processing `.cube` and `.dx` scalarfield files

- ### `edd_cube`: Calculates the electron density difference between two `.cube` files.

  Fortran program by Max Phipps.

## QM/MM

- ### `tinktep`: Interface between ONETEP and TINKER.

  A set of bash scripts by Jacek Dziedzic. Allows for performing QM/MM calculations with
  ONETEP as the QM engine and TINKER as the MM engine. Supports GAFF and the AMOEBA
  polarisable force-field. Mutual polarisation between the QM and MM subsystems is possible.
  An example is included.

## Other

- ### `prep_optados_eels`: Creates the files necessary to execute optados after a conduction run.

  Bash script by Nicholas Hine.

## Dependencies

Other scripts included here may depend on these being in your `PATH`. Put them there.

- ### `linear_fit`: Calculates a least-squares fit of f(x)=ax+b to XY, 2-column data.

  Returns *a*, *b*, their standard errors, and the regression coefficient.

- ### `getcol`: Returns the *n*-th column of a file.
- ### `getline`: Returns the *n*-th line of a file.
- ### `findavg`: Finds the average value of a column of numbers.
- ### `findmax`: Finds the maximum value of a column of numbers.
- ### `findmin`: Finds the minimum value of a column of numbers.
- ### `findrmsavg`: Finds the root-mean-square average of a column of numbers.
- ### `findsstddev`: Finds the standard deviation of a column of numbers.
- ### `findabsmax`: Finds the maximum absolute value of a column of numbers.
- ### `format_as`: Prints its argument in a given printf-like format.

