# libraries needed : 

library(vcfR)        # extraire les informations stockées dans un VCF
library(tidyverse)   # collection de packages R pour la manipulation et la visualisation des données,
library(reshape2)    # transform les data frame
library(ggplot2)     # visualisation
library(UpSetR)      # upset plot -> graphique d'intersection
library(venn)        # diagramme de venn

setwd("/Users/gaellebustarret/Projects/Stage_3A/QTL_analysis")
vcf_file <- "variant_calling/MIL2-KAD1_ParentCalling_DP_OnlyPolymorphic_NoIndel.vcf.gz"

vcf <- vcfR::read.vcfR(file = vcf_file)

#help(slotNames)
slotNames(vcf)

# the input of "vcfR2tidy" is a vcfR object
# to know the object type : is(vcf)

# to keep only specific section : 
# vcf.tidy.list <- vcfR2tidy(vcf.snv, format_fields = c("GT", "AD", "DP"), info_fields = c("AC", "AN", "MQ")) 

# you can download a csv file to use for QTLseqr here

vcf.tidy.list <- vcfR2tidy(vcf)

names(vcf.tidy.list)

head(vcf.tidy.list$meta)

head(vcf.tidy.list$fix)

head(vcf.tidy.list$gt)

# "ChromKey", "POS" common columns
vcf.tidy <- dplyr::full_join(vcf.tidy.list$fix , vcf.tidy.list$gt, by = c("ChromKey", "POS"))
head(vcf.tidy)

## exmple :
# CHROM <- "Chr1"
# POS <- 10000
# REF <- "A"
# ALT <- "T"
# ID  <- paste(CHROM, POS, REF, ALT, sep="_") # résultats : Chr1_10000_A_T

# For each variant, if ID is not null, then creat one like the example
vcf.tidy <- dplyr::mutate(vcf.tidy, ID = dplyr::if_else(is.na(ID), paste(CHROM, POS, REF, ALT, sep="_"), ID))
head(vcf.tidy)

#stringr::str_detect(string = vcf.tidy$ALT, pattern = ",")
#dplyr::filter(vcf.tidy, stringr::str_detect(string = vcf.tidy$ALT, pattern = ",")) %>% select(ID, ALT, gt_AD) %>% head

#looking for variant with more than 1 ALT allele (varaint are separated by ',' looking for this parrtern in the column)
vcf.tidy.multiAlt <- dplyr::filter(vcf.tidy, str_detect(string = ALT, pattern = ","))
vcf.tidy.multiAlt %>% dplyr::select(ID, REF, ALT, gt_AD) %>% head()

# keep variants with only one alternative allele 
filter(vcf.tidy, !str_detect(string = ALT, pattern = ",")) %>% dplyr::select(ID, ALT, gt_AD) %>% head()

vcf.No.multiAlt <- dplyr::filter(vcf.tidy, !str_detect(string = ALT, pattern = ","))

# separating the column in two :  
vcf.No.multiAlt <- tidyr::separate(data = vcf.No.multiAlt, col = gt_AD, c("gt_ref.AD", "gt_alt.AD"), sep = ",", remove=FALSE)

#convert to numeric format
vcf.No.multiAlt <- vcf.No.multiAlt %>% dplyr::mutate(gt_ref.AD = as.numeric(gt_ref.AD), gt_alt.AD = as.numeric(gt_alt.AD))

dplyr::select(vcf.No.multiAlt, ID, Indiv, gt_AD, gt_ref.AD, gt_alt.AD) %>% head()

#vcf.No.multiAlt <- dplyr::mutate(vcf.No.multiAlt, gt_AR = gt_alt.AD/gt_DP)
#vcf.No.multiAlt <- dplyr::mutate(vcf.No.multiAlt, gt_AR = if_else(gt_DP == 0, 0, gt_alt.AD/gt_DP))

#If the raw read depth (gt_DP) is =0 then the ratio is also =0 (no division by 0)  
vcf.No.multiAlt <- dplyr::mutate(vcf.No.multiAlt, gt_AR = if_else(gt_DP == 0, 0, gt_alt.AD/gt_DP))

dplyr::select(vcf.No.multiAlt, ID, Indiv, gt_AD, gt_ref.AD, gt_alt.AD, gt_AR) %>% head()

p <- ggplot(data=vcf.No.multiAlt) + 
  geom_histogram(aes(QUAL, fill=Indiv), bins = 35) + 
  xlab("Quality score") + xlim(c(0,100)) +  theme_minimal()
p

# gt_DP : Number of high-quality bases
p <- ggplot(data=vcf.No.multiAlt) + 
  geom_histogram(aes(gt_DP, fill=Indiv), bins = 35) + 
  xlab("Variant coverage")+ xlim(c(0,300))  + theme_minimal()
p

p <- ggplot(data=vcf.No.multiAlt) + 
  geom_histogram(aes(gt_DP, fill=Indiv), bins = 35) + 
  xlab("Variant coverage")+ xlim(c(0,100))  + theme_minimal()
p

p <- ggplot(data=vcf.No.multiAlt) + 
  geom_histogram(aes(gt_alt.AD, fill=Indiv), bins = 35) +
  xlab("Alternative allele depth") + xlim(c(0,50)) + theme_minimal()
p

p <- ggplot(data=vcf.No.multiAlt) + 
  geom_histogram(aes(gt_AR, fill=Indiv), bins = 35) + 
  xlab("Alternative allele ratio") + theme_minimal()
p

# regarder plus en detail l'indiv outlier
#table(vcf.No.multiAlt$Indiv)
#table(vcf.No.multiAlt$gt_DP,vcf.No.multiAlt$Indiv)

# Removed variants with quality score below 30
vcf.No.multiAlt.Q30 <- dplyr::filter(vcf.No.multiAlt, QUAL >= 30)

# Position covered by less than 4 reads are removed
# Alternative allele need to be coverred by more than 2 reads

vcf.No.multiAlt.flt <- dplyr::mutate(vcf.No.multiAlt.Q30, gt_DP = if_else(gt_DP < 4, 0, as.numeric(gt_DP)),
                                     gt_alt.AD = if_else(gt_alt.AD < 2, 0, as.numeric(gt_alt.AD)))

p <- ggplot(data=vcf.No.multiAlt.flt) + 
  geom_histogram(aes(QUAL, fill=Indiv), bins = 35) + 
  xlab("Quality score") + xlim(c(0,100)) +  theme_minimal()
p

p <- ggplot(data=vcf.No.multiAlt.flt) + 
  geom_histogram(aes(gt_DP, fill=Indiv), bins = 35) + 
  xlab("Variant coverage")+ xlim(c(0,100))  + theme_minimal()
p

p <- ggplot(data=vcf.No.multiAlt.flt) + 
  geom_histogram(aes(gt_alt.AD, fill=Indiv), bins = 35) +
  xlab("Alternative allele depth") + xlim(c(0,50)) + theme_minimal()
p

p <- ggplot(data=vcf.No.multiAlt.flt) + 
  geom_histogram(aes(gt_AR, fill=Indiv), bins = 35) + 
  xlab("Alternative allele ratio") + theme_minimal()
p

## Pour calculer une fraction de count entre plusieurs valeurs sur les diagrammes
# Calculate counts for the specified range
count_range <- sum(vcf.No.multiAlt.flt$gt_AR >= 0.3 & vcf.No.multiAlt.flt$gt_AR <= 0.75)
 
# Calculate total count
total_count <- length(vcf.No.multiAlt.flt$gt_AR)
 
# Calculate fraction
fraction <- count_range / total_count
fraction

# on peut aussi appliquer des filtre sur l'allèle ratio

# Transform the data matrix to "wide" -> ID in row et Indiv in column
vcf.mat.AR <- reshape2::dcast(vcf.No.multiAlt.flt, ID~Indiv, value.var="gt_AR")
head(vcf.mat.AR)

# Remove the ID column
vcf.mat.AR <- dplyr::select(vcf.mat.AR, -ID)

# Replace NA by 0
vcf.mat.AR <- dplyr::mutate_at(vcf.mat.AR, .vars = names(vcf.mat.AR),
                               .funs = function(x){dplyr::if_else(is.na(x), 0, as.double(x))})
head(vcf.mat.AR)

# Venn function needs a binary matrix (absence / presence)
vcf.mat.bin <- dplyr::mutate_at(vcf.mat.AR, .vars = names(vcf.mat.AR), .funs = function(x){dplyr::if_else(x > 0, 1, 0)})

#head(vcf.mat.bin)

venn::venn(vcf.mat.bin, zcolor = "style")

UpSetR::upset(vcf.mat.bin, sets=names(vcf.mat.AR),
              order.by = "freq", nintersects = NA,
              mainbar.y.label = "SNP Intersections", 
              sets.x.label = "Number of SNP")

UpSetR::upset(vcf.mat.bin, sets=names(vcf.mat.AR),
              keep.order = T, nintersects = NA,
              mainbar.y.label = "SNP Intersections", 
              sets.x.label = "Number of SNP")

# ## liste correspondance espece-genotype
# 
# # nom des genotypes
# table(vcf.tidy.list$gt$Indiv)
# 
# # vecteur des especes 
# Mdomestica <- c("X0342", "X1301", "X2953", "X6918", "X7759", "X8201", "X8396", "X8710", "X8737", "X9267")
# Msieversii <- c("MAL0944", "MAL1062")
# Msylvestris <- c("MAL0084", "MAL0900")
# 
# vcf.mat.AR <- reshape2::dcast(vcf.tidy.flt, ID~Indiv, value.var="gt_AR")
# # remplacer les NA par des 0 
# vcf.mat.AR <- dplyr::mutate_at(vcf.mat.AR, .vars = names(vcf.mat.AR),
#                                            .funs = function(x){dplyr::if_else(is.na(x), 0, as.double(x))})
# head(vcf.mat.AR)
# vcf.mat.bin <- dplyr::mutate_at(vcf.mat.AR, .vars = names(vcf.mat.AR), 
#                                             .funs = function(x){dplyr::if_else(x > 0, 1, 0)})
# 
# 
# # ajout de 3 colonnes a la matrice vcf.mat.bin
# vcf.mat.bin <-  dplyr::mutate(vcf.mat.bin, M_domestica = rowSums(dplyr::select(vcf.mat.bin,one_of(Mdomestica))), 
#                               M_sieversii = rowSums(dplyr::select(vcf.mat.bin,one_of(Msieversii))), 
#                               M_sylvestris = rowSums(dplyr::select(vcf.mat.bin,one_of(Msylvestris))))
# head(vcf.mat.bin)
# 
# 
# scp.mat.bin <-dplyr::filter(vcf.mat.bin,M_domestica ==length(Mdomestica), M_sieversii==length(Msieversii), M_sylvestris==length(Msylvestris))
# 
# head(scp.mat.bin)

# UpSetR::upset(scp.mat.bin, sets=c(Mdomestica,Msieversii,Msylvestris),
#               order.by = "freq", nintersects = NA,
#               mainbar.y.label = "SNV Intersections", 
#               sets.x.label = "Number of SNV")
# 
# UpSetR::upset(scp.mat.bin, sets=names(scp.mat.bin),
#               keep.order = T, nintersects = NA,
#               mainbar.y.label = "SNV Intersections", 
#               sets.x.label = "Number of SNV")

## Creating a table for QTLseqr

library(dplyr)
library(tidyr)

# Group by CHROM, POS, REF, and ALT
grouped <- vcf.tidy %>%
     group_by(CHROM, POS, REF, ALT)

# Select necessary columns
grouped <- select(grouped, CHROM, POS, REF, ALT, Indiv, gt_PL, gt_DP, gt_AD)

# Pivot the data to wide format
result <- grouped %>%
    pivot_wider(names_from = Indiv, 
                values_from = c("gt_PL", "gt_DP", "gt_AD"),
                names_sep = "_")

# Rename the columns
result <- result %>%
  rename("MIL2-KAD1_E.PL" = "gt_PL_MIL2-KAD1_E_BamQualFilt-f3Q30F264_nodup.bam") %>%
  rename("MIL2-KAD1_L.PL" = "gt_PL_MIL2-KAD1_L_BamQualFilt-f3Q30F264_nodup.bam") %>%
  rename("MIL2-KAD1_E.DP" = "gt_DP_MIL2-KAD1_E_BamQualFilt-f3Q30F264_nodup.bam") %>%
  rename("MIL2-KAD1_L.DP" = "gt_DP_MIL2-KAD1_L_BamQualFilt-f3Q30F264_nodup.bam") %>%
  rename("MIL2-KAD1_E.AD" = "gt_AD_MIL2-KAD1_E_BamQualFilt-f3Q30F264_nodup.bam") %>%
  rename("MIL2-KAD1_L.AD" = "gt_AD_MIL2-KAD1_L_BamQualFilt-f3Q30F264_nodup.bam")

# Writing a csv file
write.csv(result, "variant_calling/MIL2-KAD1_ParentCalling_DP_OnlyPolymorphic_NoIndel.csv", row.names = FALSE)

# Writing .table file
write.table(result, "variant_calling/MIL2-KAD1_ParentCalling_DP_OnlyPolymorphic_NoIndel.table", sep = "\t", row.names = FALSE)