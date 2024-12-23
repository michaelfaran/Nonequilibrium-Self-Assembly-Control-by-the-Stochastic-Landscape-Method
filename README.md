# Nonequilibrium Self-Assembly Control by the Stochastic Landscape Method

This repository contains the post-processed data and MATLAB scripts used to reproduce all figures presented in the manuscript and Supplementary Information (SI) of:

**Title**: *Nonequilibrium Self-Assembly Control by the Stochastic Landscape Method*  
**Authors**: Michael Faran and Prof. Gili Bisker

---

## Purpose
This repository supports the reproducibility of results and visualizations from the study. It includes:
### **Post-Processed Data**
This repository contains post-processed data aggregated and organized by figures, as follows:

- **Energy and Distance Trajectories**:
  - Examples of energy and distance versus time trajectories (Fig1).
  
- **Time to First Assembly Distributions (Tfas)**:
  - Data for Tfas distributions (Fig2, Fig1SI, Fig2SI).

- **Simulation Ensemble Results**:
  - Aggregated energy trajectory segments versus time and corresponding measures, used to construct the Stochastic Landscape (Fig3D, Fig4SI, Fig5SI).

- **Energy Trajectories and Trends**:
  - Example energy trajectories and their trends (Fig3BCEF, Fig3SI).

- **Assembly Yield and Control Protocol Analysis**:
  - Data showing the assembly yield and Tfas under different simulation conditions with the control protocol (Fig4, Fig7SI, Fig8SI, Fig9SI).

- **Equilibrium and non Equilibrium Comparisons**:
  - Tfas histograms with and without equilibrium examples (Fig6SI).

- **Drive Activation Analysis**:
  - Average distance before, during, and after drive activation under different conditions (Fig5).

- **Visual Data**:
  - Additional visual data, such as snapshots of individual simulations, used to create figures.

This organization ensures all data necessary to reproduce the figures in the manuscript and Supplementary Information is readily accessible.

- **MATLAB scripts**: `.m` files to process data and recreate the manuscript's figures.

This repository extends the work on the **Stochastic Landscape Method (SLM)**, which has been previously applied and implemented in the following studies:
1. [**Non-equilibrium Self-Assembly Time Forecasting by the Stochastic Landscape Method**](https://pubs.acs.org/doi/full/10.1021/acs.jpcb.3c01376): Repository available at [SA_UI](https://github.com/OkTAU16/SA_UI).
2. [**A Stochastic Landscape Approach for Protein Folding State Classification**](https://pubs.acs.org/doi/full/10.1021/acs.jctc.4c00464): Repository available at [Protein Folding Classification](https://github.com/michaelfaran/A-Stochastic-Landscape-Approach-for-Protein-Folding-State-Classification).

---

## Repository Structure
- **Main Folder**:
  - MATLAB scripts (`.m` files) for generating figures.
  - Each script corresponds to a specific figure, e.g., `Fig2.m` generates **Figure 2** in the main text, while `Fig2_SI.m` generates **Figure 2** in the SI.

- **Data Subfolders**:
  - Each figure script has an associated `data` subfolder containing all necessary processed data files.

---

## How to Use
1. Clone the repository to your local machine.
2. Open MATLAB.
3. Run the desired script from the main folder to generate its corresponding figure.
   - Example: Running `Fig2.m` will reproduce **Figure 2** in the main text.
4. Ensure the `data` subfolders are correctly placed alongside their corresponding scripts.

---

## Requirements
- **MATLAB**: Version R2020b or later is recommended.
- All required data files are provided in the `data` subfolders.

---

## Contact
For any questions or clarifications, please contact:  
**Michael Faran**  
michaelfaran@gmail.com
faranmic@mail.tau.ac.il

