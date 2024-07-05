# to install a new package, run install.packages("qtl")

#load the package
library(QTLseqr)
library(ggplot2)
library(dplyr)
library(tidyr)

#Set sample and file names
HighBulk <- "MIL2-KAD1_L"
LowBulk <- "MIL2-KAD1_E"
setwd("/Users/gaellebustarret/Projects/Stage_3A/QTL_analysis")
table_file <- "variant_calling/MIL2-KAD1_ParentCalling_DP_OnlyPolymorphic_NoIndel.table"
csv_file <- "variant_calling/MIL2-KAD1_ParentCalling_DP_OnlyPolymorphic_NoIndel.csv"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- paste0(rep("", 5), 1:5)

# Chroms <- paste0(rep("NC_00307", 5), c("0.9", "1.7", "4.8", "5.7", "6.8"))

#Import SNP data from file
df <- importFromGATK(
  file = table_file,
  highBulk = HighBulk,
  lowBulk = LowBulk,
  chromList = Chroms
)

#Filter SNPs based on some criteria
df_filt <- filterSNPs(
  SNPset = df,
  refAlleleFreq = 0.20,
  minTotalDepth = 100,
  maxTotalDepth = 400,
  minSampleDepth = 40,
  minGQ = 99
)

#Run G' analysis
df_filt <- runGprimeAnalysis(SNPset = df_filt, 
                             windowSize = 1e6,
                             outlierFilter = "deltaSNP")

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis( 
      SNPset = df_filt, 
      windowSize = 1e6, 
      popStruc = "F2",
      bulkSize = 250, 
      replications = 10000, 
      intervals = c(95, 99)
)

# Convert data to long format for ggplot2
df_long <- df_filt %>%
  pivot_longer(cols = starts_with("SNPindex"), 
               names_to = "phenotype", 
               values_to = "SNP_index") %>%
  mutate(phenotype = ifelse(phenotype == "SNPindex.HIGH", "Late", "Early")) %>%
  mutate(POS_Mb = POS / 1000000)
                             
# Create the plot
p <- ggplot(df_long, aes(x = POS_Mb, y = SNP_index, color = phenotype)) +
  geom_smooth(method = "auto") +
  facet_wrap(~ CHROM, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("Late" = "skyblue", "Early" = "orange")) +
  labs(x = "Position (Mb)", y = "SNP Index", color = "phenotype") +
  theme(legend.position = c(.9, .85))

# Data frame containing marker positions, chromosomes, and names
markers <- data.frame(
  CHROM = c("4", "4", "4", "4", "4", "4", "5", "5", "5", "5"),
  POS_Mb = c(0.2702020, 1.1259825, 1.2412465, 5.7254595, 9.0143655, 9.1975435, 3.1764150, 5.1716235, 5.3453005, 5.8284525),
  name = c("FRI", "LD", "GA1", "CRY1", "ESD4", "PHYD", "FLC", "CO", "FRL1", "TFL2"),
  line_number = 1:10
)

# Add vertical lines and labels to the plot
p + 
  geom_vline(data = markers, aes(xintercept = POS_Mb, linetype = as.factor(line_number)), color = "black") +
  geom_label(data = markers, aes(x = POS_Mb, y = max(df_long$SNP_index), vjust = c(0, 1.1, 2.2, 0, 0, 1.1, 0, 1.1, 2.2, 3.3), label = line_number),
             label.padding = unit(0.1, "lines"), size = 2.5, 
             fill = "white", color = "black") +
  scale_linetype_manual(values = rep("dashed", 10), labels = paste(1:10, markers$name, sep = ": "), name = "Markers") +
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  facet_wrap(~ CHROM, scales = "free_x", nrow = 1) +
  theme(legend.position = "right")

# Build the plot to extract the data
plot_data <- ggplot_build(p)

# Extract the smoothed data
smooth_data <- plot_data$data[[1]]

# Specify the x value and chromosome you are interested in
x_value <- 8.3573750 # replace with your specific x value
chromosome <- "5" # replace with your specific chromosome value

# Filter the smoothed data for the specific x value and chromosome
smoothed_values <- smooth_data %>%
  filter(PANEL == chromosome) %>%
  group_by(group) %>%
  filter(abs(x - x_value) == min(abs(x - x_value))) %>%
  ungroup()

# Display the smoothed values
print(smoothed_values)

#Plot
plotQTLStats(
  SNPset = df_filt,
  var = "deltaSNP", 
  plotIntervals = TRUE)

plotQTLStats(
  SNPset = df_filt,
  var = "Gprime", 
  plotThreshold = TRUE, 
  q = 0.01
)

plotQTLStats(
  SNPset = df_filt,
  var = "negLog10Pval", 
  plotThreshold = TRUE, 
  q = 0.01
)

setwd("/Users/gaellebustarret/Projects/Stage_3A/QTL_analysis/R_analysis")

#export summary CSV
getQTLTable(
  SNPset = df_filt,
  alpha = 0.01,
  export = TRUE,
  fileName = "MIL2-KAD1_ParentCalling_F2_QTLs.csv"
)