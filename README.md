# Elmasri_GRIN2B
Data and code for:

Elmasri M., Aziz W., Karachaliou E., Steele O.G., Sakimura, K. and Penn, A.C. Synaptic dysfunction by mutations in GRIN2B and GRIN2A: Subunit identity dominating over mutation classification. *In preparation*

## Table of contents

### Custom MATLAB functions used to analyse recording traces

**./mdocs**
- combiRec.m: Analyses AMPA-EPSCs and NMDA-EPSC responses from an I-V experiment  
- NMDARec.m: Analyses NMDA-EPSCs measured +20 mV (in prescence of NBQX)  
- ephysIO.m: Loads electrophysiology files  
- wcp.m: Analysis of test pulse to calculate whole-cell recording properties  
- rscomp.m: Off-line series resistance compensation  
- chebexp.m: Chebyshev algorithm for fitting exponential decays  

### Response data used for statistical analysis
Units for peak, decay, charge, rise, dt50, fwhm, GluN1 and Homer1c are pA, ms, pC, ms, ms, ms, au and au respectively 

**./data**
- n2b_dko_mutant_nmdar.dat
- n2b_mutant.dat
- n2b_ko_nmdar.dat  
- n2ab_C436R.dat

### R markdown (and knitted HTML) documents containing code and output relating to statistical analysis of the above data

**./rdocs**
- n2b_dko_mutant_nmdar_peak.Rmd  

- n2b_dko_mutant_nmdar_charge.Rmd  
 
- n2b_dko_mutant_nmdar_decay.Rmd
 
- n2b_mutant.Rmd  

- n2ab_C436R.Rmd 

- n2b_ko_nmdar.Rmd   

  
### Vector graphics files of graphs created in R for manuscript figures
 
./img
