


# Proteomics Workflow

## Data



## Scripts

1. `01-data_import.Rmd` - Imports and log2 transforms SomaScan protein matrix counts in addition to metadata.

2. `02-protein_lookup_table.Rmd` - Uses SomaScan.db to map Protein PROBEIDs to gene symbol and Uniprot ID. Compresses mean PROBEID value for gene symbols with >1 PROBEID. Outputs protein matrix with gene symbols and generates a protein lookup table to map PROBEID to gene symbol.

3. `03-patient_info_table.Rmd` - Subsets protein matrix counts and metadata to same 105 patients as transcriptomics workflow. Defines workflow parameters

4. `04-fold_assignment.Rmd` - Uses rsample package to split samples into 5 folds for cross-validation. Guarantees a minimum of 6 No Evidence samples per each fold to ensure balanced representation across folds.

5. `05-lasso_cv.Rmd` - Uses LASSO logistic regression to select proteins per each train/test split. Options selecting lambda which yields 1 feature per fold, 3 features per fold, or cross-validation LASSO at lambda.1se (cross_validation error within 1 standard error of lambda.min).

6. `06-random_forest.Rmd` - Selected proteins are trained in a random forest model (ntrees = 10000). For each split, out-of-fold LRTI probabilities are generated for held-out test samples.

7. `07-roc_curve_auc.Rmd` - Receiver Operating Characteristic (ROC) curve is generated from the out-of-fold lasso/RF probabilites of Definite LRTI status. AUC is reported per-fold with DeLong 95% Confidence Interval.

8. `08-differential_abundance.Rmd` - Uses limma to generate Differential Abundance (DA) volcano plot from protein matrix counts. Generates data frame from DA results which includes significance values and log2FC based on LRTI adjudication, with sex and age as covariates.

9. `09-pathway_analysis_gsea.Rmd` - Performs gene set enrichment analysis (GSEA) on ranked Differential Abundance results to identify biologically relevant pathways. Enrichment is computed using ClusterProfiler package and hallmark and reactome gene sets are pulled from MSigDB. Running enrichment score plots and dot plots are generated for select top hallmark and reactome pathways. 


































# Transcriptomic LRTI-classifier

## Scripts

1. `01-data_import.Rmd` - Imports metadata file and subsets to same 105 patients as proteomic workflow. Defines workflow parameters

2. `02-fold_assignment.Rmd` - Uses rsample package to split samples into 5 folds for cross-validation. Guarantees a minimum of 6 No Evidence samples per each fold to ensure balanced representation across folds.

3. `03-lasso_vst.Rmd` - Uses LASSO logistic regression to select genes per each train/test split. Options selecting lambda which yields 1 feature per fold, 3 features per fold, or cross-validation LASSO at lambda.1se (cross_validation error within 1 standard error of lambda.min). VST (variance stabilzed transform) parameters are estimated using only training samples, while test fold is transformed separately using frozen training parameters.

4. `04-random_forest.Rmd` - Selected genes are trained in a random forest model (ntrees = 10000). For each split, out-of-fold LRTI probabilities are generated for held-out test samples.

5. `05-roc_curve_auc.Rmd` - Receiver Operating Characteristic (ROC) curve is generated from the out-of-fold lasso/RF probabilites of Definite LRTI status. AUC is reported per-fold with DeLong 95% Confidence Interval.

6. `06-differential_expr.Rmd` - Uses limma + voom method to generate Differential Expression (DE) volcano plot from host gene counts. Generates data frame from DE results which includes significance values and log2FC based on LRTI adjudication, with sex and age as covariates.

7. `07-rna_protein_conc_plot.Rmd` - Imports DE results from transcriptomic workflow and DA results from proteomics workflow to generate an RNA/Protein concordance plot, mapping log2FC of proteins vs. genes based on LRTI adjudication. Linear model fit to jointly significant genes/proteins.

8. `08-rna_protein_loglog_plot.Rmd` - Maps protein vs. gene abundance from protein matrix (log2 transformed) and vst gene matrix.

9. `09-pathway_analysis_gsea.Rmd` - Performs gene set enrichment analysis (GSEA) on ranked Differential Expression results to identify biologically relevant pathways. Enrichment is computed using ClusterProfiler package and hallmark and reactome gene sets are pulled from MSigDB. Running enrichment score plots and dot plots are generated for select top hallmark and reactome pathways. Differential Abundance results are imported from proteomics workflow to generate a dot plot comparing top pathways in both protein and gene sets.

10. `10-supplemental_gsea_plots.Rmd` - Generates supplementary running enrichment score plots for additional top pathways.

11. `11-viral_rpm_IFNs.Rmd` - Plots standardized abundance of top IFNs (selected using GSEA) for both protein and RNA against viral RPM. Protein abundance is imported from protein matrix (log2 transformed) and vst gene matrix. Viral RPM data is obtained from `microbe_reports_bgfilter.csv` and subsetted to same 105 patients used in both workflows. Viral RPM is computed by sum viral RPM from all reported viruses. Additionally plots z-scored composite IFN score against viral RPM. Composite IFN score is calculated from summing scaled mean expression for both protein and RNA counts from IFN gene sets determined by GSEA.

12. `12-biomarker_concordance_plot.Rmd` - For select biomarkers, plots z-scored protein matrix (log2 tansformed) counts against z-scored vst gene matrix counts, where each dot is one patient_id (n=105). Additional plots color by LRTI adjudication and fit a linear model to each LRTI adjudication.


