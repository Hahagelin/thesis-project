## Overview of code

- qc.Rmd: Quality control of data.
- abspqn.Rmd: Normalization of NPX data using the "norm.AbsPQN.R" script.
- evaluation_abspqn.Rmd: Evaluates and checks the difference between the two different versions of the AbsPQN algorithm.
- table_one.Rmd: Creates "Table 1" from the thesis, an overview of the demographics of the dataset.
- visualization.Rmd: Creates PCA and UMAP plots, among other things to get an overview of the data.
- association_unadjusted.Rmd: Association analysis between protein levels and clinical variables, also adjusting for them. Also contains: Levene's test, other test testing for regional differences.
- association_adjusted.Rmd: Association analysis of protein and clinical variables after adjusting for age and BMI.
- feature_engineering.Rmd: Creates protein-protein ratio data frames.
- model_hub.Rmd: Trains, tests, and tunes models on NPX data to predict risk of breast cancer.

- norm.AbsPQN.R: Script which contains the algorithm for AbsPQN.
- association_functions.R: Contains the association analysis functions.
- model_framework.R: Contains the different model algorithms used.
