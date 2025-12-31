# Packages
library(OrgMassSpecR)
library(foreach)
library(pmartR)
library(seqinr)
library(tidyverse)

# Set diectory and load helper functions
here::i_am("Code/generate_synthetic_peptides.R")

# Functions to simulate mass spec proteomics
digest_proteins = function(sequence, enzyme = "trypsin") Digest(sequence, enzyme = enzyme)$peptide

# given a sequence of proteins, generates an overall abundance for each
generate_protein_abundances = function(){
  
}

# given a sequence of proteins, generates a coefficient for peptide j in protein i
generate_peptide_effects = function(){
  
}

# given a sequence of proteins and a group vector, generates a coefficient for protein i in group k
generate_group_effects = function(){
  
}

# provided a sequence of proteins and a group vector, combines the protein, peptide, and group effects to form the e_data
combine_effects = function(){
  
}


sequence = "MATVVVEATEPEPSGSIANPAASTSPSLSHRFLDSKFYLLVVVGEIVTEEHLRRAIGNIELGIRSWDTNLIECNLDQELKLFVSRHSARFSPEVPGQKILHHRSDVLETVVLINPSDEAVSTEVRLMITDAARHKLLVLTGQCFENTGELILQSGSFSFQNFIEIFTDQEIGELLSTTHPANKASLTLFCPEEGDWKNSNLDRHNLQDFINIKLNSASILPEMEGLSEFTEYLSESVEVPSPFDILEPPTSGGFLKLSKPCCYIFPGGRGDSALFAVNGFNMLINGGSERKSCFWKLIRHLDRVDSILLTHIGDDNLPGINSMLQRKIAELEEEQSQGSTTNSDWMKNLISPDLGVVFLNVPENLKNPEPNIKMKRSIEEACFTLQYLNKLSMKPEPLFRSVGNTIDPVILFQKMGVGKLEMYVLNPVKSSKEMQYFMQQWTGTNKDKAEFILPNGQEVDLPISYLTSVSSLIVWHPANPAEKIIRVLFPGNSTQYNILEGLEKLKHLDFLKQPLATQKDLTGQVPTPVVKQTKLKQRADSRESLKPAAKPLPSKSVRKESKEETPEVTKVNHVEKPPKVESKEKVMVKKDKPIKTETKPSVTEKEVPSKEEPSPVKAEVAEKQATDVKPKAAKEKTVKKETKVKPEDKKEEKEKPKKEVAKKEDKTPIKKEEKPKKEEVKKEVKKEIKKEEKKEPKKEVKKETPPKEVKKEVKKEEKKEVKKEEKEPKKEIKKLPKDAKKSSTPLSEAKKPAALKPKVPKKEESVKKDSVAAGKPKEKGKIKVIKKEGKAAEAVAAAVGTGATTAAVMAAAGIAAIGPAKELEAERSLMSSPEDLTKDFEELKAEEVDVTKDIKPQLELIEDEEKLKETEPVEAYVIQKEREVTKGPAESPDEGITTTEGEGECEQTPEELEPVEKQGVDDIEKFEDEGAGFEESSETGDYEEKAETEEAEEPEEDGEEHVCVSASKHSPTEDEESAKAEADAYIREKRESVASGDDRAEEDMDEAIEKGEAEQSEEEADEEDKAEDAREEEYEPEKMEAEDYVMAVVDKAAEAGGAEEQYGFLTTPTKQLGAQSPGREPASSIHDETLPGGSESEATASDEENREDQPEEFTATSGYTQSTIEISSEPTPMDEMSTPRDVMSDETNNEETESPSQEFVNITKYESSLYSQEYSKPADVTPLNGFSEGSKTDATDGKDYNASASTISPPSSMEEDKFSRSALRDAYCSEVKASTTLDIKDSISAVSSEKVSPSKSPSLSPSPPSPLEKTPLGERSVNFSLTPNEIKVSAEAEVAPVSPEVTQEVVEEHCASPEDKTLEVVSPSQSVTGSAGHTPYYQSPTDEKSSHLPTEVIEKPPAVPVSFEFSDAKDENERASVSPMDEPVPDSESPIEKVLSPLRSPPLIGSESAYESFLSADDKASGRGAESPFEEKSGKQGSPDQVSPVSEMTSTSLYQDKQEGKSTDFAPIKEDFGQEKKTDDVEAMSSQPALALDERKLGDVSPTQIDVSQFGSFKEDTKMSISEGTVSDKSATPVDEGVAEDTYSHMEGVASVSTASVATSSFPEPTTDDVSPSLHAEVGSPHSTEVDDSLSVSVVQTPTTFQETEMSPSKEECPRPMSISPPDFSPKTAKSRTPVQDHRSEQSSMSIEFGQESPEQSLAMDFSRQSPDHPTVGAGVLHITENGPTEVDYSPSDMQDSSLSHKIPPMEEPSYTQDNDLSELISVSQVEASPSTSSAHTPSQIASPLQEDTLSDVAPPRDMSLYASLTSEKVQSLEGEKLSPKSDISPLTPRESSPLYSPTFSDSTSAVKEKTATCHSSSSPPIDAASAEPYGFRASVLFDTMQHHLALNRDLSTPGLEKDSGGKTPGDFSYAYQKPEETTRSPDEEDYDYESYEKTTRTSDVGGYYYEKIERTTKSPSDSGYSYETIGKTTKTPEDGDYSYEIIEKTTRTPEEGGYSYDISEKTTSPPEVSGYSYEKTERSRRLLDDISNGYDDSEDGGHTLGDPSYSYETTEKITSFPESEGYSYETSTKTTRTPDTSTYCYETAEKITRTPQASTYSYETSDLCYTAEKKSPSEARQDVDLCLVSSCEYKHPKTELSPSFINPNPLEWFASEEPTEESEKPLTQSGGAPPPPGGKQQGRQCDETPPTSVSESAPSQTDSDVPPETEECPSITADANIDSEDESETIPTDKTVTYKHMDPPPAPVQDRSPSPRHPDVSMVDPEALAIEQNLGKALKKDLKEKTKTKKPGTKTKSSSPVKKSDGKSKPLAASPKPAGLKESSDKVSRVASPKKKESVEKAAKPTTTPEVKAARGEEKDKETKNAANASASKSAKTATAGPGTTKTTKSSAVPPGLPVYLDLCYIPNHSNSKNVDVEFFKRVRSSYYVVSGNDPAAEEPSRAVLDALLEGKAQWGSNMQVTLIPTHDSEVMREWYQETHEKQQDLNIMVLASSSTVVMQDESFPACKIEL"

xx = Digest(sequence)


sequence1 = OrgMassSpecR::example.sequence # Human Serum Albumin
unique_peps = unique(digest_proteins(sequence1))
log_scale_peptide_effects = rnorm(length(unique_peps), sd = 0.1)
pepdat = data.frame(Peptide = unique_peps, Value = 1000 * exp(log_scale_peptide_effects))
digest_n_proteins(sequence1, 500, log_scale_peptide_effects)




# Question with this linear model framework: A PTM will alter some peptides sequences, so some of the simuated coefficients for peptides are no longer valid. How to handle this?
    # Simulate coefficients for the new sequences?





##






































# Pick some sequences
sequence1 = OrgMassSpecR::example.sequence # Human Serum Albumin
sequence2 = "MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"
sequences = c(sequence1, sequence2)

G1_mean_abundances = c(1000, 2000)
G2_mean_abundances = c(1000, 1000)

expected_FC = data.frame(Protein = sequences, Expected_FC = log2(G1_mean_abundances / G2_mean_abundances))

X1 = generate_n_samples(sequences = sequences, 
                        mean_abundances = G1_mean_abundances, 
                        n_samples = 2, 
                        noise_fn = function(peptide_counts) peptide_counts, 
                        samp_id_base = NULL, 
                        enzyme = "trypsin")


X2 = generate_n_samples(sequences = sequences, 
                        mean_abundances = G2_mean_abundances, 
                        n_samples = 2, 
                        noise_fn = function(peptide_counts) peptide_counts, 
                        samp_id_base = NULL, 
                        enzyme = "trypsin")

X = full_join(X1, X2, by = 'Peptide')

# Generate peptide data
pep_e_data = X
pep_f_data = data.frame(Samp = names(pep_e_data)[names(pep_e_data) != "Peptide"], group = rep(c("1", "2"), each = 2))
pep_e_meta = get_peptide_mapping(sequences)

# Make into pmart object
pepdat = pmartR::as.pepData(e_data = pep_e_data,
                            f_data = pep_f_data,
                            e_meta = pep_e_meta,
                            edata_cname = "Peptide",
                            emeta_cname = "Protein",
                            fdata_cname = "Samp")

# Standard transformation and median centering
pepdat = pmartR::edata_transform(pepdat, "log2")
pepdat = pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = F)
pepdat = pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
pepdat = pmartR::group_designation(pepdat, main_effects = 'group')

# Rollup
methods = c('rollup', 'rrollup', 'zrollup')
combine_fns = c("mean", "median")
combinations = expand.grid(method = methods, combine_fn = combine_fns)
out_dat = c()
for(i in seq_len(nrow(combinations))){
  method = combinations$method[i]
  combine_fn = combinations$combine_fn[i]
  single_pep = method == "zrollup"
  
  protdat = pmartR::protein_quant(pepdat, method, combine_fn = combine_fn, single_pep = single_pep)
  tmp_res = imd_anova(protdat, test_method = 'anova')
  res = data.frame(tmp_res %>% 
                     select(Protein, starts_with("Fold_change")) %>%
                     left_join(., expected_FC, by = 'Protein'), 
                   Method = method, 
                   Combine_Function = combine_fn)
  
  out_dat = rbind(out_dat, res)
}
out_dat = as_tibble(out_dat)

ggdat = out_dat %>% mutate(SEL = (Fold_change_1_vs_2 - Expected_FC)^2)
ggplot(ggdat, aes(x = Method, y = SEL, color = Combine_Function)) + 
  geom_boxplot() + 
  ylab("Squared Error Loss of log2FC Recovery") + 
  labs(color = 'Combine Function') + 
  theme_bw()