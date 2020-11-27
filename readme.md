This Python module quantifies peak areas and assign them to a species
based on known retention times and finds mole composition provided a
sensitivity factor for FID or TCD chromatograms.

This is a FOSS alternative to the ChemStation and other softwares
provided by Agilent. Provided the encoder/decoder, it may be extended to
chromatographs of any type and from any manufacturer. It may be possible
to avoid entirely the use of proprietary software by using the serial
communications I/O on, for example, the 6890, but these hardware and
data acquisition issues are not here discussed.

# Data I/O

A decoder and encoder for Agilent .CH files is provided. 

Please do not try to use the encoder to fake data records. It only
changes the bits at a few positions in the file (it might not
even be readable by Agilent software--I developed this outside any
environment with Agilent ChemStation software) and there are likely
checksums in the header to ensure the data reported at the end of the
file is not corrupted. 

## Sensitivity Factors

There are simple ways to calculate sensitivity factors for flame
ionization based on the functional groups in a molecule. This requires
parsing structural formula in computer formats such as InChI or
Canonical Smiles, as atom attributes such as "carbonyl carbon" are not
stored in databases like PubChem. Such parsers are available in all
languages including python (see PySmiles, or chemplexity's molecule
repository), and the definition of these functional groups does not
require looking beyond the second nearest neighbors. For example, a
carbonyl group being any oxygen which is double bonded to an oxygen, the
acyl groups and carboxyl groups both being supersets of the carbonyl
group which require second nearest neighbor information.

You can even use pubchem RESTful API to obtain an abstract syntax tree
like representation (nested dictionaries) of the chemical structure, or
in a set of data file formats they support, from any chemical structure
specifier they support (cloud-based parsing).

See "Technical Note: Development of chemoinformatic tools to enumerate
functional groups in molecules for organic aerosol characterization"
by Giulia Ruggeri and Satoshi Takahama for a set of chemoinformatic
definitions of functional groups in terms of SMARTS patterns. In fact
the RDKit C++ (with python API) package allows enumeration of functional
groups for molecules specified with SMILES, SMARTS, or other ASCII
structural formula standards. It is merely a matter of obtaining the
counts of heteroatoms (N and O in addition to C) and functional groups
to evaluate sensitivity factors. See "Calculation of flame ionization
detector relative response factors using the effective carbon number
concept". For thermal conductivity detectors it is required to calculate
the gas phas thermal conductivity. In fact many models exist for doing
this, some of which are enumerated even in the foundational text from the 1960s
"Transport Phenomena" by Bird, Stewart, and Lightfoot, and it is one
of the cases where theory works exceptionally well because the thermal
conductivity is in the ideal gas limit only determined by pairwise
molecular interactions in their collision cross-sections. TCDs are
necessarily low in sensitivity and used for simple molecules in large
amounts, such as hydrogen and water vapor for which thermal conductivity
data is readily available.

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
a peak is sufficient to subtract the baseline. Some smoothing routines
are provided in literature and supported here, but they can introduce
artifacts such as negative signal.

## Peak Finding and Area Integration

Peak finding is a subset of signal processing and is a mature
field. Scipy offers a convenient API which merely requires articulating
what you want found, without specifying how to find it. For integration,
again due to the high time resolution, the simplest integration method
of trapezoidal rule is as accurate as higher order methods such as
Simpson's rule up to integrating Legendre polynomials fits in intervals
using the Chebyshev nodes (Gauss-Legendre quadrature).

## Peak Deconvolution

In addition to a wide literature on the subject, there are at
least two highly accessible blog posts, one from Chris Ostrouchov
(`https://chrisostrouchov.com/post/peak_fit_xrd_python/`) and another
from Emily Graceripka (`http://www.emilygraceripka.com/blog/16`). The
problem is a minimization of the cumulative square error by the fitting
parameters for a set of profiles (Gaussian, Lorentzian, and so on)
which the user initializes. An algorithm for testing which number of
profiles and of which type, with a penalty on more profiles, may also
be developed. Since most chromatograms analyzed in laboratories are
repetitive coming from a DOE on the same chemistry, it may suffice to
use a software such as `fityk` for obtaining a good guess to be used for
all chromatograms. Using lmfit a simple example is given here in the
`deconv.py` script. Arguably if peaks are convoluted in chromatography,
the solution is to change the elution temperature profile and gas
flow rates, unlike in some measurements where peak convolution is
physically necessary. However there are cases where the chemical
identities are sufficiently similar and the required data collection
interval sufficiently long that, out of practical considerations, peak
deconvolution is required.

# Codes Used and References Cited

- The binary file format parser for Agilent GCs from the `chemplexity` MATLAB package
- Whittaker smoothing algorithm (from pip package `whittaker_smoother`)
- Asymmetric least squares baseline subtraction algorithm (from literature)
- The SciPy peak finding algorithms 
- The curve fitting functions (available freely elsewhere) from the `rampy` python package
