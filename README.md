# Multi-Metric Comparison of Machine Learning Imputation Methods with Application to Breast Cancer Survival

## About

This repository contains all the necessary code and resources to reproduce the results of our study on machine learning imputation methods, using a breast cancer survival dataset as a working example. We evaluated several single and multiple imputation methods to handle missing data, crucial in clinical prognostic studies. Our comprehensive assessment covered various performance metrics each with different perspective, including Gower’s distance, estimation bias, and predictive accuracy, among others.

## Simulation Experiment

### Dataset

- **Simulated Dataset**: Created with 30% Missing At Random (MAR) values to mimic real-world scenarios of missing data in clinical studies.
  
### Imputation Methods Evaluated

- **Single Imputation (SI) Methods**: KNN, missMDA, CART, missForest, missRanger, missCforest.

- **Multiple Imputation (MI) Methods**: miceCART, miceRF.

### Performance Metrics
- Gower’s distance, estimation bias, empirical standard error, coverage rate, length of confidence interval, predictive accuracy, proportion of falsely classified (PFC), normalized root mean squared error (NRMSE), AUC, and C-index scores.

## Repository Structure

```
.
└── R/
    ├── utils.R                    # Utility functions for data generation, imputation, and performance evaluation
    └── simulation_script.R        # Main script to run the simulations and imputation experiments
```


## Citation

If you use the code in this repository in your research, please cite the repository as the manuscript is currently **under consideration**.
