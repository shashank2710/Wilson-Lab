PUPRISA
PUmp-PRobe Image Stack Analysis
https://gitlab.com/wilsonjwcsu/puprisa

Getting started: run puprisa/puprisa.m



Moved to gitlab on 4/27/2015.



Edit files in the MATLAB directory, then use the batch file to keep this copy up-to-date.

Standards directory:
This directory contains pump/probe spectroscopic standards for eu/pheo imaging. This should be expanded some day to contain other pigments.

Command-line utilities:

* puprisa_concatenateStacks.m
  
  This function concatenates (joins) several delay stacks together. This is useful
  if a delay stack needs to be constructed of more time delay points than the 
  image scanner software will support. 

* puprisa_mosaic.m

  Constructs a mosaic .dat file from a number of image stacks acquired at different
  x and y coordinates.
