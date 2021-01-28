library(bbbq)

haplotype <- "HLA-A-01:01"
n <- 9
protein_sequence <- "MNWKVLEHVPLLLYILAAKTLILCLTFAGVKMYQRKRLEAKQQKLEAERKKQSEKKDN"
topology <- "0000000000111111111111111111111000000000000000000000000000"
t <- tibble::tibble(
  epitope = create_n_mers(protein_sequence, n),
  topology = create_n_mers(topology, n)
)
t
t$ic50 <- 0.0
t$ic50 <- bbbq::predict_ic50s(
  protein_sequence = protein_sequence,
  peptide_length = n,
  haplotype = haplotype,
  ic50_prediction_tool = "EpitopePrediction"
)$ic50
epiprepreds::get_ic50_threshold(
  peptide_length = n,
  haplotype_name = haplotype_name,
  percentile = 0.02
)
t$ic50 <- round(t$ic50)
t
readr::write_csv(t, "~/smimm11a.csv")
