# Reproduction Materials for "A Model Selection Approach to Interference"

This repository contains materials to reproduce results in the titled paper.

See [README.pdf](https://github.com/gustavo-diaz/model_interference_rep/blob/main/README.pdf) for a copy of these instructions.

Visit https://gustavodiaz.org/files/research/dle_test.pdf for paper.

Visit https://gustavodiaz.org/files/research/dle_test_appendix.pdf for appendix.

# Reproduction
Source `makefile.R` to reproduce the analyses. Outputs pdf files of figures and tables.

# Manifest

- `makefile.R`: Reproduces analyses
- `application.R`: Creates figures and table for reanalysis of experiment in Ghana analyzed in section 4
- `make_sims.R`: Simulates experiments as described in section 5, sourcing this file is optional and will take days on a personal computer
- `visualize_sims.R`: Creates figures for section 5
- `GHAadjmatrix.csv`: Road adjacency network of electoral areas in Ghana. Used in `make_sims.R`
- `IchinoSchuendeln_spilloversJOP.RData`: Ghana experiment used in `application.R`
- `IchinoSchuendeln_spilloversJOP_ELAdistances.RData`: Distance matrix of ELAs in Ghana. Used in `application.R`
- `sim_df.rda`: Output from `make_sims.R`. Used in `visualize_sims.R`

# Environment

## Hardware

### For `application.R` and `visualize_sims.R`:

OS: Windows 10 Education 64-bit (10.0, Build 19043)  
Processor name: AMD Ryzen 3 2200G    
Processor speed: 3.5GHz (4 CPUs)  
Memory: 16 GB

### For `make_sims.R`:

OS: GNU/Linux distribution Scientific Linux 6.3  
Nodes: 1  
CPUs per node: 24  
Memory per CUP: 8 GB

This is the compute cluster at the School of Earth, Society, and Environment at the University of Illinois. See <http://www.econ.uiuc.edu/~lab/the-cluster.html> for details.

## Software

```
R version 4.0.5 (2021-03-31) -- "Shake and Throw"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)
```

## Packages

`future.apply_1.6.0`  
`future_1.18.0`  
`caret_6.0-86`   
`lattice_0.20-41`  
`interference_0.1.0` (https://github.com/szonszein/interference)  
`sna_2.5`  
`network_1.16.0`  
`statnet.common_4.3.0`  
`GGally_2.1.0`  
`reshape2_1.4.4`  
`MESS_0.5.6`  
`sparsereg_1.2`  
`MASS_7.3-53`  
`glmnet_3.0-2`  
`Matrix_1.2-18`  
`DeclareDesign_0.28.0`  
`estimatr_0.30.2`  
`fabricatr_0.14.0`  
`randomizr_0.20.0`  
`broom_0.7.2`  
`data.table_1.14.2`  
`forcats_0.5.0`  
`stringr_1.4.0`  
`dplyr_1.0.7`  
`purrr_0.3.4`  
`readr_1.4.0`  
`tidyr_1.1.4`  
`tibble_3.1.5`  
`ggplot2_3.3.5`  
`tidyverse_1.3.0`  
`here_0.1`  

## Run time (hh:mm:ss)

`application.R`: 00:01:01

`make_sims.R`: 11:11:14

`visualize_sims.R`: 00:00:10
