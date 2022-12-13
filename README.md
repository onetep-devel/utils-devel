# Utilities
Utility scripts and other resources for the ONETEP linear-scaling DFT code.
Download this straight into the `utils` directory of your ONETEP installation.

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

  *Note*: This script needs another script to be present in your `PATH`: `summarise`.

- ### `kern_vs_ngwf`: Helps investigate convergence indicators of a ONETEP single point energy run.

  Bash/awk script by Nicholas Hine. 
  Views the convergence of an ongoing ONETEP calculation, and displays various coloured
  indicators as to the convergence of the calculation. Examine the results to view
  decreases in energy, commutator, and ratio of optimisation energy gain from kernel and 
  NGWF optimisation and spot common problems (requires some experience to use successfully).

- ### `occ_check`: Diagnoses any potential issues with occupancies.

  Bask/awk script by Jacek Dziedzic. Reads a ONETEP `.out` file and scans the final bandstructure
  output near the Fermi level to see if it follows the Aufbau principle, essentially reporting
  on any spurious electrons above the band gap. This is meant to be used only for insulating
  systems -- in metals we obviously expect non-fractional occupancies.

  *Note*: This script needs another script to be present in your `PATH`: `getcol`.

## Format conversion

- ### `dat2xsf`: Converts `.dat` files for ONETEP into `.xsf` files for direct loading into xcrysden.

  Bash/awk script by Nicholas Hine. See the header of this script for detailed instructions.

- ### `geom2xsf`: Converts `.geom` or `.md` files produced by ONETEP into `.xsf` files for loading into xcrysden.

  Bash/awk script by Nicholas Hine. See the header of this script for detailed instructions.

- ### `geom2xyz`: Converts `.geom` or `.md` files from ONETEP or CASTEP into `.xyz` for animation in xmol.

  Perl script by Greg Pearce. This is from 2002 and no longer maintained. It might not function correctly in this day and age.

- ### `dat2xyz`: Converts the coordinates in a ONETEP `.dat` file to a `.xyz` file.

  See below, under 'Operations on input (`.dat`) files'.

- ### `out2charge`: Extracts atomic charges (Mulliken, NPA/NBO or DDEC) from a ONETEP `.out` file. 

  Bash/awk script by Jacek Dziedzic. Extracts atomic charges (Mulliken, NPA/NBO or DDEC) from a ONETEP `.out` file. 
  Additionally, if an `.xyz` file exists, it is used, together with the charges to generate a `.vtf` file 
  that is readily understood by VMD and where 'Coloring method -> charge' simply works. This can
  be used to easily visualise any kinds of atomic charges in VMD.

  Bash/awk script by Jacek Dziedzic. Extracts atomic charges (Mulliken, NPA/NBO or DDEC) from a ONETEP `.out` file. 

- ### `pdb2dat`: Creates a ONETEP input (`.dat`) file from a `.pdb` file.

  Bash/awk script by Jacek Dziedzic. Converts the coordinates and atom types from a `.pdb`
  file to a ONETEP input (`.dat`) file. ONETEP keywords to go into the `.dat` file are 
  read from separate file -- `pdb2dat.header`, so that they can be changed without having 
  to modify the script itself. If `pdb2dat.header` is absent, the script looks for
  it in the directory the script itself is located in. Box size and atom translation 
  can be adjusted with the command line arguments. Run this script with no arguments for help.

- ### `xyz2dat`: Creates a ONETEP input (`.dat`) file from an `.xyz` file.

  Bash/awk script by Jacek Dziedzic. Converts the coordinates and atom types from a `.xyz`
  file to a ONETEP input (`.dat`) file. ONETEP keywords to go into the `.dat` file are 
  read from separate file -- `xyz2dat.header`, so that they can be changed without having 
  to modify the script itself. If `xyz2dat.header` is absent, the script looks for
  it in the directory the script itself is located in. Box size and atom translation 
  can be adjusted with the command line arguments. Run this script with no arguments for help.

- ### `usp2upf`: Converts a pseudopotential in the CASTEP `.usp` form to unified pseudopotential format.

  Fortran program by PWSCF group. Distributed under the terms of the GNU General Public license.

## Molecular dynamics

- ### `md2thermo`: Calculates thermodynamical quantities from a ONETEP `.md` file.

  Bash/awk script by Jacek Dziedzic. Can be used to process an MD trajectory (`.md` file) to get statistical information on energy conservation and temperature.

  *Note*: This script needs several other scripts to be present in your `PATH`. These are: `getcol`, `getline`, `findmin`, `findmax`, `findavg`, `findrmsavg`, `findabsmax`, `findstddev`, `format_as`, and `linear_fit`. These scripts are also provided here.

- ### `md2vel`: Extracts velocities from a ONETEP `.md` file.

  Bash/awk script by Jacek Dziedzic. Processes a `.md` trajectory file and extracts
  atomic velocities into a form compatible with a multi-frame `.xyz` file.

- ### `md_clean_duplicates`: Cleans up an `.md` trajectory from multiple restarts.

  Bash/awk script by Jacek Dziedzic. Reads an `.md` trajectory file, copying input to output, except for 
  duplicate steps (those with identical times), which are removed. Only the first instance of any time is retained. 
  This is useful for processing `.md` files combined from multiple runs, where global history has been used, 
  leading to short time intervals present at the end of one run, and then again at the beginning of a restared run.

## Manipulation of input (`.dat`) files

- ### `dat2bounds`: Shows information about the geometry of a ONETEP `.dat` file.

  Bask/awk script by Jacek Dziedzic. Shows information about the geometry of a ONETEP `.dat` file,
  including margins to box sides and recommended settings for cutoff Coulomb. This is a simple,
  limited script -- e.g. it only works for orthorhombic boxes, it expects the box to be specified 
  in atomic units, not angstroem; it also assumes 7a0 NGWFs, but it's easy to manually correct the results.

- ### `dat2xyz`: Converts the coordinates in a ONETEP `.dat` file to an `.xyz` file.

  Bask/awk script by Jacek Dziedzic. Reads the `%positions_abs` and `%lattice_cart` blocks,
  and produces an `.xyz` file corresponding to the system defined by the `.dat` file. This is
  a simple, limited script -- e.g. it only works for orthorhombic boxes, it expects the box to be specified 
  in atomic units, not angstroem.

  *Note*: This script needs another script to be present in your `PATH`: `getcol`.

- ### `datcentre`: Centres the system in a `.dat` file.

  Bash/awk script by Jacek Dziedzic. Translates all atoms in a `.dat` file to move the
  geometric centre of the system to the centre of the box.  This is
  a simple, limited script -- e.g. it only works for orthorhombic boxes, it expects the box to be specified 
  in atomic units, not angstroem.

  *Note*: This script needs another script to be present in your `PATH`: `dat2bounds`.

- ### `datshift`: Moves the system in a `.dat` file by an offset.

  Bash/awk script by Jacek Dziedzic. Translates all atoms in a `.dat` file by offsets
  specified separately for all Cartesian directions. Leaves the box unchanged.
  This is a simple, limited script -- e.g. it only works for orthorhombic boxes, it expects the 
  box to be specified in atomic units, not angstroem.

## Processing `.cube` and `.dx` scalarfield files

- ### `edd_cube`: Calculates the electron density difference between two `.cube` files.

  Fortran program by Max Phipps. Similar functionality, and more, is offered by `dxsum`.

- ### `dx2cube`: Converts a `.dx` file to a `.cube` file.

  C++ program by Jacek Dziedzic. Binary available on request from J.Dziedzic@soton.ac.uk. 
  Converts data in a `.dx` file to a `.cube` file. Optionally reads a `.dat` file to extract also atomic positions
  to be included in the `.cube` file. Handles requisite unit conversions and sign
  convention changes automatically.

- ### `dxaxpy`: Computes *a*\**x*+*y* for data *x* in a `.dx` file, outputting results to a `.dx` file.

  C++ program by Jacek Dziedzic. Binary available on request from J.Dziedzic@soton.ac.uk. 
  Provides an easy way to scale and shift the values in `.dx` files.

- ### `dxchop`: Zeroes values above or below a threshold in a `.dx` file.

  C++ program by Jacek Dziedzic. Binary available on request from J.Dziedzic@soton.ac.uk.

- ### `dxcoarsen`: Coarsens (keeps only every $n$-th point) in a `.dx` file.

  Bash/awk script by Jacek Dziedzic. Provides an easy way of making `.dx` files
  more lightweight should they prove to be too large to handle, e.g. by VMD. Typically
  used to reduce the resolution by a factor of 2 along every direction.

- ### `dxgeom`: Provides details about the geometry of a `.dx` file.

  Bash/awk script by Jacek Dziedzic.

- ### `dxstats`: Provides statistical information about data in a `.dx` file.

  Bash/awk script by Jacek Dziedzic. Prints the number of datapoints, min and max values,
  average, sum, sum of absolute values, average absolute value, sum of squares and
  the rms value. Useful for detecting anomalies in `.dx` files.

- ### `dxintegrate`: Integrates the quantity in a `.dx` file.

  C++ program by Jacek Dziedzic. Binary available on request from J.Dziedzic@soton.ac.uk.

- ### `dxproduct`: Multiplies data in two `.dx` files, producing a `.dx` file as output.

  C++ program by Jacek Dziedzic. Binary available on request from J.Dziedzic@soton.ac.uk. 
  Computes y=*a*\**x1*\**x2*+*b* for data *x1* in one `.dx` file and data *x2* in another `.dx` file,
  with *a* and *b* being constants specified as arguments.

- ### `dxsum`: Adds data in two `.dx` files, producing a `.dx` file as output.

  C++ program by Jacek Dziedzic. Binary available on request from J.Dziedzic@soton.ac.uk. 
  Computes y=*a1*\**x1+*a2*\**x2* for data *x1* in one `.dx` file and data *x2* in another `.dx` file,
  with *a1* and *a2* being constants specified as arguments. Note how with *a1*==1 and *a2*==-1
  you can calculate differences.

- ### `dxsection`: Generates a 1D cross-section along one of the Cartesian directions of `.dx` file.

  Bash/awk script by Jacek Dziedzic. Generates a 1D cross-section along *X*, *Y* or *Z*
  from a `.dx` file. Very useful for producing plots.

  *Note*: This script needs `dxserver` to be present in your `PATH`.

## DFTB

- ### `dftb`: Parameter files for the DFTB method.

## QM/MM

- ### `tinktep`: Interface between ONETEP and TINKER.

  A set of bash scripts by Jacek Dziedzic. Allows for performing QM/MM calculations with
  ONETEP as the QM engine and TINKER as the MM engine. Supports GAFF and the AMOEBA
  polarisable force-field. Mutual polarisation between the QM and MM subsystems is possible.
  An example is included.

## Other

- ### `JTH_v1.0_PAW_datasets`: PAW dataset library produced by the developers of ABINIT, converted to a format suitable for ONETEP.

  `.tar.gz` archives containing `.abinit` and `.atompaw` files.

- ### `DDEC_atomic_densities`: Reference ion densities for DDEC3 (density derived electrostatic and chemical) charge analysis.

  `.coreconf` and `.refconf` files by Thomas Manz and Nidia Gabaldon Limas.

- ### `vdW_params`: Van der Waals parameters and a program for their determination for ONETEP's vdwcorrection module.

  Fortran source and data files by Quintin Hill.

- ### `prep_optados_eels`: Creates the files necessary to execute optados after a conduction run.

  Bash script by Nicholas Hine.

- ### `unitconv`: Converts units.

  Bash/awk script by Jacek Dziedzic. Helps with unit conversion. Examples of use:

  *What is the conversion factor from bohr to angstroem?*

  ``unitconv distance_au distance_A``

  *How many rydbergs are there in a kcal/mol?*

  ``unitconv energy_kcal_per_mol energy_ry``

  *How much is 2.5 debye in SI units (C*m)?*

  ``unitconv dipole_debye dipole_si :2.5``

  *Convert forces in `forces.txt` from a0/Ha to eV/A:*

  ``unitconv force_au force_eV_over_A forces.txt``

  For a list of supported quantities and units, run the script with no parameters.

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

