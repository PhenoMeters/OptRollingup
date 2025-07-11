here::i_am("Code/protein_sequence_simulation.R")

source(here::here("Code", "ptm_utils.R"))

lip_data = readr::read_tsv(here::here("Data", "double_pept_lytic_sites.tsv"))

all_peptides = unlist(lapply(lip_data$Sequence, \(x) strsplit(x, split = '; ')[[1]]))

unique_peptides = unique(all_peptides)

# Lets look at the frequency of each amino acid in the real data
amino_acids = unlist(lapply(unique_peptides, \(x) strsplit(x, split = "")[[1]]))
barplot(table(amino_acids), xlab = "Amino Acids", ylab = "Frequency", main = "Amino Acid Frequency Plot")

amino_acid_distribution = create_amino_acid_distribution(amino_acids)
synthetic_protein = generate_protein(1000, amino_acid_distribution)

barplot(table(synthetic_protein), xlab = "(Synthetic) Amino Acids", ylab = "Frequency", main = "(Synthetic) Amino Acid Frequency Plot")
