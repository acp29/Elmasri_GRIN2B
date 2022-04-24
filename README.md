# Elmasri_GRIN2A
Data and code for:

Elmasri M, Hunter DW, Winchester G, Bates EE, Aziz W, Van Der Does DM, Karachaliou E, Sakimura K, Penn AC.  (2022) Common synaptic phenotypes arising from diverse mutations in the human NMDA receptor subunit GluN2A. *Commun Biol*. 5(1):174. doi: [10.1038/s42003-022-03115-3](https://www.nature.com/articles/s42003-022-03115-3) 

## Table of contents

### Custom MATLAB functions used to analyse recording traces

**./mdocs**
- combiRec_IV.m: Analyses AMPA-EPSCs and NMDA-EPSC responses from an I-V experiment  
- combiRec.m: Analyses NMDA-EPSCs measured +20 mV (in prescence of NBQX)  
- ephysIO.m: Loads electrophysiology files  
- wcp.m: Analysis of test pulse to calculate whole-cell recording properties  
- rscomp.m: Off-line series resistance compensation  
- chebexp.m: Chebyshev algorithm for fitting exponential decays  

### Response data used for statistical analysis
Units for peak, decay, charge, rise, dt50, fwhm, GluN1 and Homer1c are pA, ms, pC, ms, ms, ms, au and au respectively 

**./data**
- n2a_dko_mutant_nmdar.dat: relates to Figures 2 and S1.1-1.2  
- n2a_mutant_imaging.dat: relates to Figures 3 and S1.3  
- n2a_mutant.dat: relates to Figure 4, S4 and S1.4-1.7  
- n2a_ko_nmdar.dat: relates to Figure 5 and S1.8-1.11  
- n2a_dose_rescue_nmdar.dat: relates to Figure S3  

### R markdown (and knitted HTML) documents containing code and output relating to statistical analysis of the above data

**./rdocs**
- n2a_dko_nmdar_peak.Rmd: relates to Figure 2b  
  Double knockout (DKO, grin2a and grin2b): NMDA-EPSC peak  
  https://rpubs.com/acp29/835110  
- n2a_dko_mutant_nmdar_peak.Rmd: relates to Figure 2c and S1.1  
  GluN2A mutant rescue in DKO: NMDA-EPSC peak  
  https://rpubs.com/acp29/835112  
- n2a_dko_mutant_nmdar_decay.Rmd: relates to Figure 2d and S1.2  
  GluN2A GOF mutant rescue in DKO: NMDA-EPSC decay  
  https://rpubs.com/acp29/835113  
- n2a_mutant_imaging.Rmd: relates to Figures 3 and S1.3  
  GluN1-SEP and Homer1c-tdTomato spine fluorescence when co-expressed with GluN2A mutants  
  https://rpubs.com/acp29/835101  
- n2a_mutant.Rmd: relates to Figure 4, S4 and S1.4-1.7   
  GluN2A mutant rescue in *grin2a* knockout: EPSC properties   
  https://rpubs.com/acp29/835104  
- n2a_ko_nmdar.Rmd: relates to Figure 5 and S1.8-1.11   
  grin2a knockout: NMDA-EPSC properties  
  https://rpubs.com/acp29/835105  
- n2a_dose_rescue_nmdar.Rmd: relates to Figure S3  
  grin2a KO dose-dependent rescue with WT GluN2A: NMDA-EPSC peak and decay  
  https://rpubs.com/acp29/835115  
  
### Vector graphics files of graphs created in R for manuscript figures
 
./img
