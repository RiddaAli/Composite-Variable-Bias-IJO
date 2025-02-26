## Composite variable bias: causal analysis of weight outcomes

  **R code:**
  - 2025-02-25 NCDS-DAG models.R: it includes all the code used for the analyses including sensitivity analyses 

  **Datasets:**
  - The individual datasets of National Child Development Study (NCDS) cohort surveys at ages 23 and 33 can be downladed from the [UK Data Service](https://ukdataservice.ac.uk/).
  - The file selected_data_ncd_age23.csv contains missing values for baseline height (BaseHt). These missing values can be manually filled using the follow-up height (FUpHt) variable from the NCDS data at age 33, assuming height remains unchanged between ages 23 and 33. The updated file is then saved as selected_data_ncd_age23_imputed.csv. This will then reproduce exactly what was done for the article, otherwise there  will be small differences due to the missing values.
