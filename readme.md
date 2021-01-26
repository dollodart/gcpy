This Python module quantifies peak areas and assign them to a species
based on known retention times and finds mole composition provided a
sensitivity factor for FID or TCD chromatograms.

This is a FOSS alternative to the ChemStation and other softwares
provided by Agilent. Provided you know the encoder/decoder, it may be
extended to chromatographs of any type and from any manufacturer. It may
be possible to avoid entirely the use of proprietary software by using
the serial communications I/O on, for example, the Agilent 6890, but
these hardware and data acquisition issues are not here discussed.

# Data I/O

A decoder and encoder for Agilent .CH files is provided. The decoder was copied to python from the MATLAB Chemplexity package by James Dillon.

Please do not try to use the encoder to fake data records. It only
changes the bits at a few positions in the file. It might not
even be readable by Agilent software--I developed this outside any
environment with Agilent ChemStation software. And there are likely
checksums in the header to ensure the data reported at the end of the
file is not corrupted. 

## Sensitivity Factors

There are simple ways to calculate relative sensitivity factors
for flame ionization based on the functional groups in a set of
molecules. This requires parsing structural formula in computer formats
such as InChI or Canonical Smiles and processing that data, as atom
attributes such as "carbonyl carbon" are not stored in databases like
PubChem. Parsers are available in all languages including python: here
pysmiles is used. The definition of these functional groups does not
require looking beyond the second nearest neighbors. For example, a
carbonyl group being any oxygen which is double bonded to an oxygen, the
acyl groups and carboxyl groups both being supersets of the carbonyl
group which require second nearest neighbor information. Ethers and
esters require only second nearest neighbor information, provided the
central node used is the oxygen. The functional group determination is
coded into rules on the element of a central atom, the elements of its
neighboring (bonded) atoms, and the elements of second neighbors (atoms
bonded to atoms which are bonded to the central atom).

For thermal conductivity detectors one needs to calculate the gas phas
thermal conductivity. Models based on the kinetic theory of gases and
other theoretical and semi-empirical approaches have been used, some of
which are enumerated in pedagagical fashion as early as the 1960s in the
foundational text "Transport Phenomena" by Bird, Stewart, and Lightfoot.
It is one of the cases where theory works exceptionally well because
the thermal conductivity is in the ideal gas limit only determined by
pairwise molecular interactions in their collision cross-sections. TCDs
are necessarily low in sensitivity, in fact the thermal conductivity
is non-unique with composition for mixtures because of the highly
non-linear mixing rules for thermal conductivities. They can be used for
simple molecules like hydrogen and helium for which thermal conductiivty
data is readily available, and in systems where there are few non-trace
gases which are chemically similar so the non-uniqueness is less of a
problem.

## Retention Times

Though it should be possible to determine from some physical models the
retention time in packed columns, the modeling is evidently sufficiently
difficult that it is easier to use pure component test samples or attach
a mass spectrometer to the end of the gas chromatograph. Also the exact
characteristics of the chromatographic columns are proprietary so the
required data, e.g., to determine adsorption enthalpy, are not known.

# Quantification

## Baseline Subtraction

The time resolution and sensitivity of modern chromatographs are so
high that for most cases, a linear interpolation between the edges of
a peak is sufficient to subtract the baseline. 

## Smoothing

Some smoothing routines provided by the literature are coded here, but
they can introduce artifacts such as negative signal and generally
require inspection. One advantage to smoothing is perhaps in identifying
the functional form of the baseline so that one may fit a first or
second order polynomial. When one does asymmetric smoothing, there
results asymmetric coefficients which effectively identify the points
belonging to peaks and those to baseline, so that one can fit the
baseline only to the baseline points. But the same can be achieved with
peak finding algorithms. So smoothing is kept here only for aesthetic
purposes.

## Peak Finding 

Peak finding is a subset of signal processing and is a mature
field. Scipy offers a convenient API which merely requires articulating
what you want found, without specifying how to find it. 

## Curve Integration

Like in the case of baseline subtraction, because of the high time
resolution the simplest integration method of trapezoidal rule is as
accurate as higher order methods such as Simpson's rule, let alone
numerical integration methods like Gauss-Legendre quadrature.

## Peak Deconvolution

In addition to a wide literature on the subject, there are at
least two highly accessible blog posts, one from Chris Ostrouchov
(`https://chrisostrouchov.com/post/peak_fit_xrd_python/`) and another
from Emily Graceripka (`http://www.emilygraceripka.com/blog/16`). The
problem is a minimization of the cumulative square error by the fitting
parameters for a set of profiles (Gaussian, Lorentzian, and so on)
which the user initializes. A program for testing which number of
profiles and of which type, with a penalty on more profiles, could
be developed. Since most chromatograms analyzed in laboratories are
repetitive coming from a DOE on the same chemistry, it should suffice to
use a software such as `fityk` on a representative sample for obtaining
a good guess to be used for all chromatograms. Using `lmfit` a simple
example is given here in the `deconv.py` script. Arguably if peaks are
convoluted in chromatography, the solution is to change the elution
temperature profile and gas flow rates, unlike in some measurements
where peak convolution is physically necessary. However there are cases
where the chemical identities are sufficiently similar and the required
data collection interval sufficiently long that, out of practical
considerations, peak deconvolution is required.

# Codes Used and References Cited

- The binary file format parser for Agilent GCs from the `chemplexity` MATLAB package
- Whittaker smoothing algorithm from pip package `whittaker_smoother`
- Asymmetric least squares baseline subtraction algorithm from "Baseline Correction with Asymmetric Least Squares
  Smoothing"
- The SciPy peak finding algorithms 
- The curve fitting functions (available freely elsewhere) from the `rampy` python package

See "Technical Note: Development of chemoinformatic tools to enumerate
functional groups in molecules for organic aerosol characterization"
by Giulia Ruggeri and Satoshi Takahama for a set of chemoinformatic
definitions of functional groups in terms of SMARTS patterns.

See "Calculation of flame ionization detector relative response factors
using the effective carbon number concept" as a reference for the
effective carbon numbers.
