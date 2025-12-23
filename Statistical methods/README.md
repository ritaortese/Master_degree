**Project Description – Health Care Demand in Australia**

This project analyzes the factors that influence the number of doctor visits in Australia using a health care dataset. The main variable of interest is doctorco, which represents the number of consultations with a doctor or specialist in the previous two weeks.

The project starts with data preprocessing and exploratory data analysis (EDA) to understand the structure of the dataset, clean variables, and reduce redundancy. Several variables were combined or removed to avoid multicollinearity and improve model interpretation. Graphical analysis and correlation analysis were used to explore relationships between variables.

**Tools and software used**
R and RStudio were used for all data analysis and modeling.

Important R packages used include:

  - dplyr for data manipulation
  
  - ggplot2, GGally, corrplot for data visualization
  
  - caret for train/test split and evaluation
  
  - mgcv for Generalized Additive Models (GAM)
  
  - pscl for zero-inflated models
  
  - pROC for ROC curves and AUC in classification tasks
  
**Statistical Models used**
  - Linear Regression models, both with and without interaction terms, as a baseline.
  
  - Poisson Regression, suitable for count data like the number of doctor visits.
  
  - Generalized Additive Models (GAMs) to capture possible non-linear relationships between predictors and the response variable.
  
  - Zero-Inflated Poisson (ZIP) and Zero-Inflated Negative Binomial (ZINB) models to handle the large number of zero observations and overdispersion in the data.
  
  - Sampling techniques (undersampling, oversampling, and a combination of both) were used to address class imbalance.
  
  - Logistic Regression was used for a binary classification task to predict whether a person had zero or at least one doctor visit.

Models were compared using RMSE, AIC, ROC curves, and confusion matrices, and the results showed that more flexible models (such as GAMs and zero-inflated models) generally provided better performance.
