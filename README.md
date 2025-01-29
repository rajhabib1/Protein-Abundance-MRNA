# Predicting Protein Abundance from mRNA data

## Introduction
The R script is designed to predict protein abundance using mRNA sequences, RNA abundance data, and translation measurements. Our method integrates data preprocessing, feature extraction, and a machine learning approach utilizing a random forest model, aiming to bridge mRNA transcription and protein translation processes for biological research.
## Data Preparation
### Packages and Libraries
The script starts by installing and loading necessary R packages, which include:
- tidyverse for data manipulation,
- Biostrings for handling biological string data,
- randomForest for the machine learning model implementation,
- BiocManager for managing Bioconductor packages.
### Data Importation
The script imports several datasets for analysis:
- mRNA_sequences containing mRNA sequences in FASTA format.
- rna_cl_train, ribo_cl_train, prot_train include RNA-seq data for RNA abundance, ribosome profiling data, and protein abundance data from mass spectrometry respectively.
- region_lengths, human_info provide data on gene region lengths and experimental metadata pertaining to different human cell lines.
## Data Normalization
A custom function, normalize_data_CPM_avg, normalizes the data to Counts Per Million (CPM) and averages it across samples from the same cell line to ensure comparability by mitigating differences in sequencing depth and sample size.
## Feature Engineering
### Gene Region Lengths
The extract_gene_region_lengths function computes the lengths of 5'UTR, CDS, and 3'UTR regions from a BED file, essential for understanding mRNA regulatory effects on translation and stability.
### GC Content Calculation
Calculating the GC content for each mRNA sequence is crucial as GC-rich regions can affect the structural stability and translational efficiency of mRNA.
### mRNA Stability Index
Stability is assessed through AU-rich elements within the mRNA sequences, affecting mRNA degradation rates. This index is derived by counting AU-rich sequences and normalizing these counts by the length of the mRNA sequence.
### Kozak Sequence Analysis (Discontinued)
Initial attempts to integrate the strength of the Kozak consensus sequence surrounding the start codon as a feature were discontinued after findings suggested minimal impact on differentiating protein abundance levels across samples.
## Data Integration
The script merges all extracted features (RNA abundance, ribosome profiling data, gene region lengths, GC content, and mRNA stability) into a single dataset, aligning each gene with its corresponding features and measured protein abundance.
## Machine Learning Model
### Random Forest Implementation
The random forest algorithm is employed for its ability to handle complex, high-dimensional datasets and capture non-linear relationships among features. The dataset is formatted and integrated for input into the machine learning model.
### Model Training and Validation
The script divides the data into training and testing sets. It then trains the model and evaluates its performance using metrics such as Mean Squared Error (MSE), Root Mean Squared Error (RMSE), and R-squared:
- MSE: 2539.403
- RMSE: 50.39249
- R-squared: 0.5965595
- These metrics indicate the model’s accuracy and the variability it can explain, with an R-squared of approximately 0.60, suggesting moderate predictive power.
Even though our model’s accuracy is only 60%, we managed to add 2 additional features to predict protein abundance, which with adding more features, we could improve the accuracy.
### Cross-Validation
Cross-validation is conducted to ensure the model’s robustness and generalizability across different data subsets.
## Conclusion and Model Predictions
The final trained model is capable of making predictions on protein abundance for new RNA and translation measurement inputs. The model's predictions, demonstrated in lines such as:
- gene_id AAAS predicted_protein_abundance 167.40706
- gene_id AACS predicted_protein_abundance 159.13811

This indicates that the model not only predicts known values closely but can also extrapolate to new data. This capability proves beneficial in contexts where protein levels are pivotal, such as disease research or drug development. The distribution across various cell lines, including detailed gene-wise predictions, underscores the model's utility in predicting protein expression variability across different biological conditions.
