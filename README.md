# Python codes to compute SC functions and compute dv/v for [Viens et al. (2018, JGR)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JB015697)

## Description:
This repository contains Python codes for [Viens et al. (2018, JGR)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JB015697).
The single station cross-correlation (SC) functions for all the MeSO-net stations between February 20, 2011 to May 11, 2011 can be downloaded at: https://doi.org/10.7910/DVN/TPDPHA

This repository contains:
* The **Codes** folder:
- Plot_CC_dv_example.py: Compute dv/v measurement from the data and plot results (the data first need to be downloaded). This code uses:
  - Stretching.py: function to compute dv/v measurement with the stretching method.
  - Additional_Functions.py:
- Get_SC_functions.py: Main function used to compute single-station cross-correlation (SC) functions.

* The **Data** folder:
 - Empty folder. SC functions can be downloaded at [here](https://doi.org/10.7910/DVN/TPDPHA).
 
* The **Figures** folder:
 - To save the output of Plot_CC_dv_example.py

## Example:
The output of the **Plot_CC_dv_example.py** code when applied to the SC functions calculated at the NS7M MeSO-net station are shown below. 

![Comparaison between the different methods](https://github.com/lviens/2018_JGR/blob/master/Figures/Fig_dv_E.NS7M.png)

