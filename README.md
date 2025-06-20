# task-based-fMRI-ICN-analysis

**Activation and Dynamic Connectivity Modulations of the Brain’s Intrinsic Networks during Visuomotor, Response Control, and Working Memory Processes**

**Elif Kurt, Ali Bayram, Tamer Demiralp**

This repository contains the data, the analysis scripts and output tables used in our study investigating intrinsic connectivity networks (ICNs) during task-based fMRI. The project focuses on both intrinsic activation and dynamic functional network connectivity (dFNC) changes across working memory, visuomotor integration, and response control conditions.

## Folder Structure

- `code/`  
  Contains all custom MATLAB scripts for the analysis of ICN activation and dFNC computation.

- `data/`  
  Contains the task design file (`GLM_design_mdes.mat`) and the time series data for each subject (`ICA_timeseries_loaded_C15.mat`).
  - `GLM_design_mdes.mat`: includes the `SPM.mat` structure representing the task design matrix.  
  - `ICA_timeseries_loaded_C15.mat`: stores concatenated ICN time series for all subjects across tasks.

- `results/`  
  Includes Excel files generated from activation and connectivity analyses. These files were used for subsequent statistical analyses (e.g., in SPSS).

## Requirements

To run the analysis, the following MATLAB-based toolboxes must be installed:

- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- [GIFT Toolbox](http://trendscenter.org/software/gift/)
- [MarsBaR Toolbox](https://marsbar.sourceforge.net/) (integrated within SPM12)

> **Important Note:**  
> Please replace the default `mars_arm_call.m` file located in the MarsBaR toolbox with the modified version provided in:
>
> ```
> code/mars_arm_call.m
> ```

## Scripts and Execution

### ICN Activation Analysis

To compute intrinsic activation of ICNs across task conditions, run the following script:

```
code/ICN_intr_activ_batch.m
```

### dFNC Analysis

To compute task-modulated dynamic functional network connectivity, run:

```
code/ICN_dfnc_batch.m
```

This script internally calls the modified GIFT dFNC function:

```
code/icatb_compute_task_dfnc_HUBAL.m
```

### Output

All results are saved in Excel format within the `results/` directory.
These output files were used for statistical analyses performed in SPSS.

---

For questions or requests regarding the data and code, please contact the corresponding author. elifkurt@istanbul.edu.tr
