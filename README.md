<div align="center">

# Rotation Measure Scatter (RMS)

_✨ Code repository for Paper ✨_  
**Frequency Dependent Polarization of Repeating Fast Radio Bursts - Implications for Their Origin**

</div>

<p align="center">
  <a href="https://github.com/SukiYume/RMS">
    <img src="https://img.shields.io/badge/RotationMeasureScatter-RMS-red" alt="release">
  </a>
</p>

## Description
  
  In this paper, we report polarization measurements of five repeating FRBs. Combining these with archival observations, we identify trends of lower polarisation at lower frequencies. We model this behaviour as multi-path **Rotation Measure Scatter (RMS)**. 

  <img src="Figure/RM_Scatter_New.png" alt="RMS" width="900px" />

  Sources with higher **RMS** have higher *RM magnitude* and *scattering timescales*, indicating a complex environment near the sources, such as a supernova remnant or a pulsar wind nebula, consistent with FRBs arising from young stellar populations.

## Structure

  Figure 2, 3, and 4 of this paper are plotted in `RM-Figure.ipynb`. `Polarization-RM.ipynb` contains Figure 1 and the measurements of FRBs' Rotation Measure. Data and Figures are stored in `CalData` and `Figure`, respectively.

  `ModulePy` contains two methods for measuring RM, **RM Synthesis** in `rm_synthesis.py`, and **RM QU Fitting** in `rm_qu_fitting.py`. The code usage method shows in `code_usage.py`.

  ```bash
  `-- RMS
      |-- CalData
      |   |-- 190417-FAST.csv
      |   |-- 201124-FAST.csv
      |   |-- 201124-GBT.csv
      |   |-- RM-190303-1.txt
      |   |-- RM-190417-7.txt
      |   |-- RM-190520-437.txt
      |   |-- RM-201124-26.txt
      |   `-- RM-437.txt
      |-- Figure
      |   `-- RM_Scatter_New.png
      |-- ModulePy
      |   |-- code_usage.py
      |   |-- extract_pulse.py
      |   |-- requirements.txt
      |   |-- rm_qu_fitting.py
      |   `-- rm_synthesis.py
      |-- RM-Figure.ipynb
      |-- Polarization-RM.ipynb
      |-- LICENSE
      `-- README.md
  ```

  Usage of `rm_qu_fitting.py` - For example:

  ```python
  from rm_qu_fitting import mcmc_fit, plot_mcmc_samp

  # Q and U - 1D numpy array of a pulse, as a function of frequency.
  result, RM, RM_error_left, RM_error_right = mcmc_fit(Q, U, freq, rm_left=-10000, rm_right=10000)
  print('RM: {:.0f} +{:.0f} -{:.0f}'.format(RM, RM_error_right, RM_error_left))
  plot_mcmc_samp(result, save=False)
  ```



## TO-DO

  Modularize the measurement methods of RM.