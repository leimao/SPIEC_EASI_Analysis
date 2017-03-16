# SPIEC EASI Analysis for Microbiome Data
# Author: Lei Mao
# Date: 3/16/2017

library(devtools)
library(SpiecEasi)
library(phyloseq)
library(Matrix)

# Microbiome file name
file <- "HMPv35_100nt.biom"
# Filter percentage
parameter_filter <- 0.2
# data might be extremely big, se careful about the memory size.
data <- import_biom(file, parseFunction = parse_taxonomy_greengenes)
# Filter out the OTUs that exist in less than parameter_filter of the samples  
data_filtered <- filter_taxa(data, function(x) sum(x > 1) > (parameter_filter * length(x)), TRUE)
# Run SPIEC EASI
data_se <- spiec.easi(data_filtered, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, icov.select.params=list(rep.num = 50))
# Exact the weights of association
data_se_weights <- summary(symBeta(getOptBeta(data_se), mode = 'maxabs'))
# Save weights
write.csv(data_se_weights, file = "weights.csv", row.names = FALSE)
# Save taxonomy
write.csv(data_filtered@tax_table, file = "tax_table.csv", row.names = FALSE)


# Save log file
file_log <- "log.txt"
cat("Microbiome Data", "\n", file = file_log, append = FALSE, sep = '')
cat(as.character(file), "\n", file = file_log, append = TRUE, sep = '')
cat("Filter Parameter", "\n", file = file_log, append = TRUE, sep = '')
cat(as.character(parameter_filter), "\n", file = file_log, append = TRUE, sep = '')
cat("Number of Samples", "\n", file = file_log, append = TRUE, sep = '')
cat(as.character(dim(data@otu_table)[2]), "\n", file = file_log, append = TRUE, sep = '')
cat("Number of OTUs Before Filtering", "\n", file = file_log, append = TRUE, sep = '')
cat(as.character(dim(data@otu_table)[1]), "\n", file = file_log, append = TRUE, sep = '')
cat("Number of OTUs After Filtering", "\n", file = file_log, append = TRUE, sep = '')
cat(as.character(dim(data_filtered@otu_table)[1]), "\n", file = file_log, append = TRUE, sep = '')

