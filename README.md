# Transcriptomic LRTI-classifier

## Scripts

1. `01-data_import.Rmd` - Imports metadata file and subsets to same 105 patients as proteomic workflow. Defines workflow parameters

2. `02-fold_assignment.Rmd` - Uses rsample package to split samples into 5 folds for cross-validation. Guarantees a minimum of 6 No Evidence samples per each fold to ensure balanced representation across folds.

3. `03-lasso_vst.Rmd` - Uses LASSO logistic regression to select genes per each train/test split. Options selecting lambda which yields 1 feature per fold, 3 features per fold, or cross-validation LASSO at lambda.1se (cross_validation error within 1 standard error of lambda.min). VST (variance stabilzed transform) parameters are estimated using only training samples, while test fold is transformed separately using frozen training parameters.

4. `04-random_forest.Rmd` - Selected genes are trained in a random forest model (ntrees = 10000). For each split, out-of-fold LRTI probabilities are generated for held-out test samples.

5. `05-roc_curve_auc.Rmd` - Receiver Operating Characteristic (ROC) curve is generated from the out-of-fold lasso/RF probabilites of Definite LRTI status. AUC is reported per-fold with DeLong 95% Confidence Interval.

6. `06-differential_expr.Rmd` - Uses limma + voom method to generate differential expression (DE) volcano plot from host gene counts. Generates data frame from DE results which includes significance values and log2FC based on LRTI adjudication, with sex and age as covariates

7. `07-rna_protein_conc_plot.Rmd` - Imports DE results from transcriptomic workflow and DA results from proteomics workflow to generate an RNA/Protein concordance plot, mapping log2FC of proteins vs. genes based on LRTI adjudication. Linear model fit to jointly significant genes/proteins.

8. `08-rna_protein_loglog_plot.Rmd` - Maps protein vs. gene abundance from protein matrix (log2 transformed) and vst gene matrix

9. 

