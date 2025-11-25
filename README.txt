This document describes the code scripts and data available to replicate analyses presented in Beach et al. (2025) "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

This Rproject file contains all datasets used in analysis, aggregated at the municipio level. The main dataset used for modeling is the "mun_ds.rds" file, which has been saved as an .rds file along with a CSV and shapefile format in the /data folder. The 'data_dictionary.csv' describes each variable in the dataset. 

The following scripts included in this Rproject can be followed in order to replicate modeling and geographic analysis described in the manuscript. 

'a_covariate_exploration.R' fits reduced models to test the association between hypothesized ecological factors and gastric cancer risk. Candidate variables for fully adjusted models are identified. 

'b_fully_adjusted_model_building.R' builds multiple regression models. Fully adjusted model covariate set is selected from assessing fit statistics (WAIC and marginal likelihood) and covariate significance at 95% level. Explores spatial structure of fully adjusted model Pearson residuals and assesses multicollinearity among the fixed effects. Supplementary figures A1, A2, A7, and A8 are created in this script. Supplementary Tables A1 and A2 are created in this script.

'c_fit_models.R' fits fully adjusted FE and spatial models for all strata. Null, bivariate, and fully adjusted models are fit for each strata to assess how posterior Theta[i], adjusted spatial relative risks RR[i], Phi, and coefficient estimates may be sensitive to inclusion of each covariate. Figure 3 created in this script. Supplementary Table C1, Table C2, and Figure A6 are created in this script.

'd_priors_sensitivity_checks.R' assesses posterior distributions under different priors for coefficients and hyperparameters. Supplementary figures A3, A4, and A5 are created in this script.

'e_residual_spatial_structure.R' assesses the spatial structure of fully adjusted models (FE and spatial). Table 3 is created in this script. 

'f_mapping_visualizations.R' creates all map figures presented in the manuscript and supplementary materials. Figures 2, 4, 5 and B1-B8 created in this script.


[Index of figures/tables created]
Fig/Table - Script

Figure 2  - F
Figure 3  - C
Figure 4  - F
Figure 5  - F
Table  3  - E

Figure A1 - B
Figure A2 - B
Figure A3 - D
Figure A4 - D
Figure A5 - D
Figure A6 - C
Figure A7 - B
Figure A8 - B
Table  A1 - B
Table  A2 - B

Figure B1 - F
Figure B2 - F
Figure B3 - F
Figure B4 - F
Figure B5 - F
Figure B6 - F
Figure B7 - F

Table  C1 - C
Table  C2 - C


