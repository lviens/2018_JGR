# Python codes to compute SC functions and dv/v measurements ([Viens et al., 2018, JGR](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JB015697))

## Description:
This repository contains Python codes for [Viens et al. (2018, JGR)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JB015697).
The single station cross-correlation (SC) functions for all the MeSO-net stations between February 20, 2011 to May 11, 2011 can be downloaded at: https://doi.org/10.7910/DVN/TPDPHA

This repository contains:
* The **Codes** folder:
  - Plot_CC_dv_example.py: Compute dv/v measurement from the data and plot results (the data first need to be downloaded). This code uses:
    - Stretching.py: function to compute dv/v measurement with the stretching method.
    - Additional_Functions.py: Additional functions to plot the data.
  - Get_SC_functions.py: Main function used to compute single-station cross-correlation (SC) functions.

* The **Data** folder:
  - Empty folder. SC functions can be downloaded at [here](https://doi.org/10.7910/DVN/TPDPHA).
 
* The **Figures** folder:
  - To save the output of **Plot_CC_dv_example.py** code.

## Example:
The output of the **Plot_CC_dv_example.py** code when applied to the SC functions calculated at the NS7M MeSO-net station is shown below. SC functions are computed over 15-min data for the Z-N and Z-E components and subsequently stacked over 6 hours. The reference waveform is the stack of the SC functions between Feb. 20 to March 10th, 2011 (e.g., before the 2011 Mw 9.0 Tohoku-Oki earthquake). The SC functions are bandpass filtered between 1 and 20 Hz and the stretching is performed on the first second of the SC functions. 
* The top subplot shows the correlation coefficient between the reference waveform and each current window. A correlation coefficient equal to 1 indicates that the two waveforms are perfectly symmetric. Note that the CC for the Z-N and Z-E SC functions are averaged following [Hobiger et al. (2014, GJI)](https://academic.oup.com/gji/article/198/1/90/604971).
* The bottom subplot shows the dv/v measurements averaged over the Z-N and Z-E SC functions following [Hobiger et al. (2014, GJI)](https://academic.oup.com/gji/article/198/1/90/604971). Dark red colors indicate a high correlation coefficient and white/yellow colors indicate low correlation coefficients. The large drop of dv/v on March 11th, 2011 is caused by the Mw 9.0 2011 Tohoku-Oki earthquake. The small drop in February is caused by a precipitation event. Note that the figure below is a little bit different from that of the paper (e.g., Figure 7) as the maximum dv/v change allowed when performing the stretching is set to 10%.

![Comparaison between the different methods](https://github.com/lviens/2018_JGR/blob/master/Figures/Fig_dv_E.NS7M.png)

