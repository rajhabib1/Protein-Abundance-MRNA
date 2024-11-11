
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings", force = TRUE)
install.packages("randomForest")
library(tidyverse)
library(Biostrings)
library(stringr)
library(randomForest)

#####Loading datasets
mRNA_sequences <- readDNAStringSet("~/BioInformatics/mRNA_sequences.fa.gz")
load("~/BioInformatics/rna_cl_train.rda")
load("~/BioInformatics/ribo_cl_train.rda")
load("~/BioInformatics/prot_train.rda")
region_lengths <- read_csv("~/BioInformatics/RegionLengths.bed")
human_info <- read_csv("~/BioInformatics/human_infor_train.csv")



##########Function to normalize data using CPM and average across samples for the same cell line
normalize_data_CPM_avg = function(file_to_normalize, human_info = experiments_cell_line) {
  #Mapping experiment_alias to cell_line
  cell_line_map = setNames(human_info$cell_line, human_info$experiment_alias)
  #Renaming the columns of file_to_normalize to the corresponding cell line
  colnames(file_to_normalize) = cell_line_map[intersect(names(file_to_normalize), names(cell_line_map))]
  #Calculate total counts per sample to use for normalization
  total_counts_per_sample = colSums(file_to_normalize, na.rm = TRUE)
  #Normalize counts to CPM (counts per million)
  cpm_normalized = sweep(file_to_normalize, 2, total_counts_per_sample, FUN = "/") * 1e6
  #If multiple columns for the same cell line, average them
  unique_cell_lines = unique(human_info$cell_line)
  cpm_averaged = sapply(unique_cell_lines, function(cell_line) {
    cell_line_data = cpm_normalized[, colnames(cpm_normalized) == cell_line, drop = FALSE]
    rowMeans(cell_line_data, na.rm = TRUE)
  })
  #Convert the matrix back to a data frame and set row and column names appropriately
  cpm_averaged_df = as.data.frame(cpm_averaged)
  rownames(cpm_averaged_df) = rownames(file_to_normalize)
  colnames(cpm_averaged_df) = unique_cell_lines
  #Returning the CPM averaged data frame
  return(cpm_averaged_df)
}
#Normalizing RNA abundance data
norm_rna_abundance_cpm_avg = normalize_data_CPM_avg(rna_cl_train, human_info)
print(norm_rna_abundance_cpm_avg)
#Normalizing translation measurements data
norm_translation_measurements_cpm_avg = normalize_data_CPM_avg(ribo_cl_train, human_info)
print(norm_translation_measurements_cpm_avg)



######Extracting lengths
region_lengths_path <- "~/BioInformatics/RegionLengths.bed"
library(tidyverse)
extract_gene_region_lengths <- function(region_lengths_path) {
  # Read the file lines from the provided path
  data_lines <- read_lines(region_lengths_path)
  # Convert to a tibble and separate out the columns
  region_lengths <- tibble(line = data_lines) %>%
    separate(line, into = c("gene_id", "start", "end", "region_type"), 
             sep = "\\t", remove = TRUE, convert = TRUE) %>%
    # Calculate the length of each region
    mutate(length = end - start) %>%
    # Summarize lengths by gene and region type
    group_by(gene_id, region_type) %>%
    summarize(total_length = sum(length), .groups = 'drop') %>%
    # Reshape the summary to have CDS, UTR3, and UTR5 as separate columns
    pivot_wider(names_from = region_type, values_from = total_length)
  # Return the wide format of the length summary
  return(region_lengths)
}
length_summary_wide <- extract_gene_region_lengths(region_lengths_path)
print(length_summary_wide)


#####Extracting data
library(Biostrings)
#GC Content Percentage
calculate_gc_content_percentage <- function(mRNA_sequences) {
  library(Biostrings)
  gc_contents <- letterFrequency(mRNA_sequences, letters = c("G", "C"), as.prob = TRUE)
  gc_content <- rowSums(gc_contents)
  gc_content_percentage <- gc_content * 100
  names(gc_content_percentage) <- names(mRNA_sequences)
  gc_content_df <- data.frame(gene_id = names(gc_content_percentage), gc_content = gc_content_percentage)
  return(gc_content_df)
}
gc_content_df <- calculate_gc_content_percentage(mRNA_sequences)

#mRNA Stability (presence and count of specific sequence motifs like AU-rich elements))
#Function to calculate mRNA stability index based on a broader AU-rich elements detection
calculate_mRNA_stability <- function(mRNA_sequences) {
  #Regular expression for a broader range of AU-rich elements
  pattern <- "[AU]{2,}"  #Detects sequences of two or more A's and U's
  #Calculating the count of AU-rich elements in each sequence
  stability_scores <- sapply(as.list(mRNA_sequences), function(seq) {
    length(gregexpr(pattern, as.character(seq), perl=TRUE)[[1]][gregexpr(pattern, as.character(seq), perl=TRUE)[[1]] > 0])
  })
  #Normalizing scores by the length of the mRNA sequence to get a stability index
  stability_index <- stability_scores / width(mRNA_sequences)
  #Creating a dataframe with gene_id and their corresponding stability index
  stability_df <- data.frame(gene_id = names(mRNA_sequences), stability_index = stability_index)
  return(stability_df)
}
stability_df <- calculate_mRNA_stability(mRNA_sequences)
stability_df

#Kozak Sequence - Tried but not a good predictor
#calculate_kozak_strength <- function(mRNA_sequences) {
  #kozak_strength_scores <- sapply(as.list(mRNA_sequences), function(seq) {
    #seq <- as.character(seq)
    #Finding all ATG positions
    #start_codon_indices <- unlist(gregexpr("ATG", seq))
    #if (length(start_codon_indices) > 0) {
      #Checking for the presence of a start codon at the beginning of the sequence
      #if (start_codon_indices[1] > 6) {  # Ensure there's enough sequence before the start codon
        #Extracting the Kozak sequence context: 6 bases before to 3 bases after the ATG
        #kozak_region <- substr(seq, start_codon_indices[1] - 6, start_codon_indices[1] + 3)
        #Optimal Kozak consensus: gccRccATGG
        #Computing how many crucial positions (-3 and +4 around ATG) match this pattern
        #score <- 0
        #if (substr(kozak_region, 4, 4) == "A" || substr(kozak_region, 4, 4) == "G") { # Position -3: R = A or G
          #score <- score + 1
        #}
        #if (substr(kozak_region, 9, 9) == "G") { # Position +4: G
          #score <- score + 1
        #}
        #return(score)
      #}
    #}
    #return(0)  # No valid ATG found or not enough context
  #})
  #kozak_strength_df <- data.frame(gene_id = names(mRNA_sequences), kozak_strength = kozak_strength_scores)
  #return(kozak_strength_df)
#}
#kozak_strength_df <- calculate_kozak_strength(mRNA_sequences)
#kozak_strength_df


######Integrating data by cell line
#Defining the cell line names and initialize a list for the merged data frames
cell_lines <- c("A549", "HeLa", "HepG2", "K562", "MCF7", "U2OS")
merged_data_frames <- list()
#Preparing the length summary data by ensuring it has a 'gene_id' column
if(!"gene_id" %in% colnames(length_summary_wide)) {
  length_summary_wide$gene_id <- row.names(length_summary_wide)
}
#Looping through each cell line, merging data from multiple sources, and merging with length summary
for (cell_line in cell_lines) {
  # Create a merged data frame for the current cell line
  merged_cell <- data.frame(
    rna_abundance = norm_rna_abundance_cpm_avg[, cell_line],
    translation_measurements = norm_translation_measurements_cpm_avg[, cell_line],
    protein_abundance = prot_train[, cell_line],
    gene_id = rownames(norm_rna_abundance_cpm_avg)
  )
  #Merging with the gene length data
  merged_data <- merge(merged_cell, length_summary_wide, by = "gene_id", all.x = TRUE)
  #Storing the merged data frame in the list with the cell line name as the key
  merged_data_frames[[cell_line]] <- merged_data
  #Printing the merged data for the cell line
  print(paste("Merged data for cell line:", cell_line))
  print(head(merged_data))
}


#Merging GC content data frame with each cell line's merged data
for (cell_line in names(merged_data_frames)) {
  #Merging the GC content with the merged data frames
  merged_data_frames[[cell_line]] <- merge(merged_data_frames[[cell_line]], gc_content_df, by = "gene_id", all.x = TRUE)
  #Printing the merged data for the cell line
  print(paste("Merged data with GC content for cell line:", cell_line))
  print(head(merged_data_frames[[cell_line]]))
}
#Integrating mRNA stability into each cell line's merged data
for (cell_line in names(merged_data_frames)) {
  merged_data_frames[[cell_line]] <- merge(merged_data_frames[[cell_line]], stability_df, by = "gene_id", all.x = TRUE)
  print(paste("Merged data with mRNA stability for cell line:", cell_line))
  print(head(merged_data_frames[[cell_line]]))
}
merged_data_frames
all_data <- do.call(rbind, merged_data_frames)

#Integrating Kozak's Consensus Sequence - Discontinued
#for (cell_line in names(merged_data_frames)) {
  #merged_data_frames[[cell_line]] <- merge(merged_data_frames[[cell_line]], kozak_strength_df, by = "gene_id", all.x = TRUE)
  #print(paste("Merged data with Kozak strength for cell line:", cell_line))
  #print(head(merged_data_frames[[cell_line]]))
#}
#all_data <- do.call(rbind, merged_data_frames)




#######Training and Testing The Model
#Splitting the data into training and testing sets
set.seed(123)
training_indices <- sample(1:nrow(all_data), 0.7 * nrow(all_data))
train_data <- all_data[training_indices, ]
summary(train_data)
summary(all_data[-training_indices,])
#Check for Missing Values in the Data
summary(train_data)
#Imputing missing values in the UTR5 column with the median
train_data$UTR5[is.na(train_data$UTR5)] <- median(train_data$UTR5, na.rm = TRUE)
#Verifying that no NAs remain
sum(is.na(train_data))
test_data <- all_data[training_indices, ]
#Retrain the Random Forest Model
rf_model <- randomForest(protein_abundance ~ rna_abundance + translation_measurements +
                           CDS + UTR5 + UTR3 + gc_content + stability_index,
                         data = train_data, ntree = 500, mtry = 4, importance = TRUE)



#Check the summary of the model to ensure it has trained correctly
summary(rf_model)

#Making predictions on the test dataset
predictions <- predict(rf_model, newdata = test_data)
view(predictions)


#Calculating MSE, RMSE, and R-squared, ensuring no NAs interfere
valid_indices <- which(!is.na(predictions) & !is.na(test_data$protein_abundance))
new_predictions <- predictions[valid_indices]
new_test_data <- test_data$protein_abundance[valid_indices]
mse <- mean((new_predictions - new_test_data)^2)
rmse <- sqrt(mse)
rsq <- cor(new_predictions, new_test_data)^2
#Printing the results
print(paste("RMSE:", rmse))
print(paste("MSE:", mse))
print(paste("R-squared:", rsq))

##############Function to predict protein abundance given RNA abundance and translation information
predict_protein_abundance <- function(rna_abundance_data, translation_measurements_data) {
  #Normalizing RNA abundance and translation measurements data
  norm_rna_abundance <- normalize_data_CPM_avg(rna_abundance_data, human_info)
  norm_translation_measurements <- normalize_data_CPM_avg(translation_measurements_data, human_info)
  #Reading the mRNA sequences file
  mRNA_sequences <- readDNAStringSet("path_to_your_mRNA_sequences_fasta_file")
  #Calculating GC content from mRNA sequences
  gc_content_df <- calculate_gc_content_percentage(mRNA_sequences)
  #Calculating mRNA Stability
  stability_df <- calculate_mRNA_stability(mRNA_sequences)
  #Calculating Kozak's Sequence
  kozak_strength_df <- calculate_kozak_strength(mRNA_sequences)
  #Extracting the lengths of the 5UTR, CDS, and 3UTR from the region lengths file
  region_lengths_wide <- extract_gene_region_lengths("path_to_your_region_lengths_file")
  #Creating datasets to run the model on
  datasets <- create_datasets(norm_rna_abundance, norm_translation_measurements, region_lengths_wide, gc_content_df, stability_df, kozak_strength_df)
  #Getting protein abundance predictions
  protein_predictions <- predict(rf_model, newdata = datasets)
  return(protein_predictions)
}




##Observations for different cell lines

# Observations for A549 cell line
a549_data <- merged_data_frames[["A549"]]
summary(a549_data)


##Cross-validation

library(caret)
set.seed(123)
folds <- createFolds(all_data$protein_abundance, k = 5)
cv_results <- lapply(folds, function(fold_index) {
  train_fold <- all_data[-fold_index, ]
  test_fold <- all_data[fold_index, ]
  
  # Train the model on the fold
  fold_model <- randomForest(protein_abundance ~ rna_abundance + translation_measurements +
                               CDS + UTR5 + UTR3 + gc_content + stability_index,
                             data = train_fold, ntree = 500, mtry = 4, importance = TRUE)
  
  # Validate the model on the test fold
  fold_predictions <- predict(fold_model, newdata = test_fold)
  
  # Calculate performance metrics
  fold_mse <- mean((fold_predictions - test_fold$protein_abundance)^2)
  fold_rmse <- sqrt(fold_mse)
  fold_rsq <- cor(fold_predictions, test_fold$protein_abundance)^2
  
  return(list(MSE = fold_mse, RMSE = fold_rmse, R_squared = fold_rsq))
})

# Calculate average performance metrics across all folds
avg_mse <- mean(sapply(cv_results, function(x) x$MSE))
avg_rmse <- mean(sapply(cv_results, function(x) x$RMSE))
avg_rsq <- mean(sapply(cv_results, function(x) x$R_squared))

print(paste("Average MSE:", avg_mse))
print(paste("Average RMSE:", avg_rmse))
print(paste("Average R-squared:", avg_rsq))

