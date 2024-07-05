import pandas as pd
import os

os.chdir("/Users/gaellebustarret/Projects/Stage_3A/QTL_analysis/R_analysis")

FT_genes = ['FLM', 'MAF','PHYA', 'CRY1', 'CRY2', 'PHYB', 'PHYD', 'PHYE', 'PFT1', 'FKF1', 'GI', 'CO', 'CDF1', 'SPY', 'LD', 'FVE', 'FLK', 'FRL1', 'FRL2', 'FES1', 'FRI', 'HUA2', 'HOS1', 'PIE1', 'ESD4', 'ELF5', 'FLC', 'VIN3', 'VIN3-L', 'TFL2', 'SVP', 'GA1', 'GA', 'GID1', 'RGL1', 'RGL2', 'SLY1', 'FPF1', 'GAI', 'RGA', 'ATMYB33', 'EBS', 'AGL24', 'SOC1', 'FT', 'FD', 'FDP', 'TSF', 'TFL1', 'ATC', 'BFT', 'MFT', 'AP1', 'CAL', 'LFY', 'AP3', 'PI']

# Read the DataFrame from the CSV file
dfIntersect = pd.read_csv("MIL2-KAD1_F2_QTLs_norownames_annotation_intersect_test.csv", sep=";", dtype=str)

# Create a mapping dictionary for chromosome names
chromosome_mapping = {
    'CP002684.1': '1',
    'CP002685.1': '2',
    'CP002686.1': '3',
    'CP002687.1': '4',
    'CP002688.1': '5'
}

# Filter for rows where 'type' is 'gene'
exon_df = dfIntersect[dfIntersect['type'] == 'gene']

# Filter rows where the column "Unnamed: 14" starts with 'Note='
filtered_df = exon_df[exon_df['Unnamed: 14'].str.startswith('Note=', na=False)]

# Extract the gene and locus tag information from the filtered DataFrame
locus_tags = filtered_df['Unnamed: 11'].str.replace('ID=gene-', '')
genes = filtered_df['Unnamed: 13'].str.replace('Name=', '')
notes = filtered_df['Unnamed: 14'].str.replace('Note=', '')

# Extract the chromosome information
chromosomes = filtered_df['chrom'].map(chromosome_mapping)

# Create a DataFrame with genes, locus tags, and chromosome names
df_genes = pd.DataFrame({
    'Gene': genes,
    'Chromosome': chromosomes,
    'Locus_Tag': locus_tags,
    'Start_Position': filtered_df['pos_start'],
    'End_Position': filtered_df['pos_end'],
    'Note': notes
})

# Select rows where "Note" contains "flower"
df_genes_flower = df_genes[(df_genes['Note'].str.contains('flower', case=False, na=False)) | (df_genes['Gene'].isin(FT_genes))]

# Export the DataFrame to a CSV file
df_genes_flower.to_csv('MIL2-KAD1_F2_QTL_FT_genes.csv', index=False)
