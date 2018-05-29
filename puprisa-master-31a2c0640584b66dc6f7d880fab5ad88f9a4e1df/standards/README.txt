Standards used for linear decomposition are located in this directory. Each .mat file contains the pump-probe response of a single chromophore.

Each .mat file here has at least 3 variables:

  name:	the name of the standard

  t:	the probe delay vector, in picoseconds

  x:	the pump-probe signal

Standards with multiple wavelength combinations have these additional
variables:

  wlpairs: an array of pairs of pump/probe wavelength combinations,
           e.g. [705, 817; 725, 817]. Each row is a combination; the 1st 
           column contains pump wavelengths, the 2nd column contains probe 
           wavelengths

  tindex:  a 1D array, a list of indices to the row in wlpairs that
           corresponds to the wavelength combination used for each element
           in t.  For example, if both wavlength combinations were sampled
           with 5 different time delays, tindex would be 
           [1 1 1 1 1 2 2 2 2 2]


============================================================================
=== Single wavelength pair standards for eu/pheo ===========================

euPheo720pump810probe.mat:
Tom's cuvette data (used in the Sci. Transl. Med. publication).

euPheo_featherHair_720pump810probe.mat:
Mary Jane's black hair / red feather data (not normalized to eu/pheo concentration).


============================================================================
=== Multi-wavelength pair standards for eu/pheo ============================

eumelaninEDTA_705_817_725_817.mat:

  EDTA-washed eumelanin, recorded at both 705nm/817nm and 725nm/817nm
  pump/probe by MJS 2013-05-29. Data compiled by JWW 2013-06-05.


eumelaninFe_705_817_725_817.mat:

  Iron-loaded eumelanin, recorded at both 705nm/817nm and 725nm/817nm
  pump/probe by MJS 2013-05-29. Data compiled by JWW 2013-06-05.


pheomelanin_705_817_725_817.mat:

  Synthetic pheomelanin, recorded at both 705nm/817nm and 725nm/817nm
  pump/probe by MJS 2013-05-29. Data compiled by JWW 2013-06-05.


hemoglobin_705_817_725_817.mat:

  Hemoglobin from red blood cells in a conjunctival melanoma, recorded at
  both 705nm/817nm and 725nm/817nm pump/probe by MJS 2013-05-29. Data 
  compiled by JWW 2013-06-05.




==== additional places to find standards =====
Hemoglobin:
NREM ROIB 40x z=-55 zDelayStack from 3/5/2012
NREM ROID 40x z=+35 zDelayStack from 3/12/2012