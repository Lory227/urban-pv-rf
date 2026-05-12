# urban-pv-rf
This repository contains selected MATLAB source code from my PhD research on urban PV potential / surface albedo / radiative forcing modelling. The code is shared as a reviewable example of the modelling workflow, software structure, and data-processing logic.

The full reproducible research archive, including large LiDAR/DSM, geospatial, and associated input/output materials, is hosted on the 4TU Research Data repository:

https://data.4tu.nl/datasets/136d3808-cb58-4276-8d5b-4a4e2cd5e38c/1

Due to the size of the original spatial datasets, this GitHub repository does not include the full input data. Instead, it provides the main scripts and workflow documentation to illustrate the implementation approach.

**Author**: Yilong Zhou  
**Date**: May 2025

This document applies to the scripts created for the project studying the impacts of urban photovoltaic (PV) integration on radiative forcing. The scripts are divided into three parts:

1. Albedo simulation before and after PV integration  
2. Annual PV energy calculation  
3. Positive and negative radiative forcing calculation  

---

##  Input Data Description

All input data are stored in the folder `INPUTS`.

### i) INPUTS folder contains:

- `albedoResults.mat`: Daily, monthly, and yearly albedo results before and after PV integration  
- `astmg173.xls`: SMARTS data for Delft
- `Brightness4.mat`: Key meteorological data for the simulation, based on Meteonorm    
- `GHI_tot.mat`: Yearly GHI at daytime hours in Delft  
- `Materials.mat`: Pre-processed spectral reflectance of materials from ASTER library 
- `PVlaminate.mat`: Measured spectral reflectance of SUNPOWER IBC solar cell  
- `PVresults.mat`: Annual roof PV yield results for Delft  
- `Reflectance.xlsx`: Measured spectral reflectance of diffuse glass and IBC PV cell  

### ii) Subfolder `DSMs` contains:

- DSM TIFFs for Delft: `C_37EN1`, `C_37EN2`, `C_37EZ1`, `C_37EZ2`  
- DSM TIFFs containing adjacent tiles: `Delft1_9T`, `Delft_9T`  

### iii) Subfolder `Shapefiles` contains:

- Ecosystem Unit Map shapefiles   
- PV polygons for two integration scenarios from the PV simulation  

---

## I. Albedo Simulation

- Run `albedoCalculator.m` in `Albedo Calculation` folder  
- Select the DSM tile and whether to apply PV material integration  
- If yes, select scenario one or two  
- The results are pre-simulated and saved in `albedoResults.mat` in the `INPUTS` folder  

*Note*: Adjust number of parallel workers for your system

---

## II. Annual PV Energy Calculation

- Run `Delft_PV_Potential.m` in `PVPotential Calculation` folder  
- The code automatically calculates annual PV output for both scenarios  
- Results are saved in `PVresults.mat` in `INPUTS` folder  

*Note*: Adjust number of parallel workers for your system

*Note*: Files in the 'Coefficients' folder are used for simplified skyline based model

---

## III. Radiative Forcing Calculation

- Run `RFcalculator.m` in the `RF Calculation` folder  
- Outputs positive and negative radiative forcing results  
- To test different carbon emission intensities, edit `kwh_co2_NL` in line 231  

---

**Citation Required**: If you are using the code from this work, please cite the corresponding publication (to be updated, currently under review).

**Acknowledgements to external data and tools**: 
- SMARTS model by NREL (https://www.nrel.gov/grid/solar-resource/smarts):
	  
	--  Gueymard, C., "Parameterized Transmittance Model for Direct Beam and Circumsolar Spectral Irradiance", Solar Energy, Vol. 71, No. 5, pp. 325-346, 2001.

  --  Gueymard, C., "SMARTS, A Simple Model of the Atmospheric Radiative Transfer of Sunshine: Algorithms and Performance Assessment". Professional Paper FSEC-PF-270-95. Florida Solar Energy Center, 1679 Clearlake Road, Cocoa, FL 32922, 1995.

- ASTER library: (doi:10.1016/j.rse.2008.11.007)  
- Meteonorm data: (https://mn8.meteonorm.com/en/meteonorm-version-8) 
- LiDAR data: (https://www.pdok.nl/introductie/-/article/actueel-hoogtebestand-nederland-ahn)
- Building footprints: (https://www.pdok.nl/introductie/-/article/basisregistratie-adressen-en-gebouwen-ba-1)
- Ecosystem Unit Map: (https://www.cbs.nl/en-gb/background/2017/12/ecosystem-unit-map)
- Radiative forcing calculation: 

	-- Wohlfahrt Georg, Tomelleri Enrico and Hammerle Albin (2021) “The albedo-climate penalty of hydropower reservoirs”. 
	   Zenodo. doi: 10.5281/zenodo.4432576.
---

For questions, contact: Yilong Zhou (musezyl@gmail.com)

