library(bbbq)

n <- 9
haplotype <- "HLA-A-01:01"
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
t$ic50 <- round(t$ic50)
ic50_threshold <- epiprepreds::get_ic50_threshold(
  peptide_length = n,
  haplotype_name = haplotype,
  percentile = 0.02
)
t$is_binder <- FALSE
t$is_binder <- t$ic50 < ic50_threshold

readr::write_csv(t, "~/smimm11a_is_binder.csv")


n_haplotypes <- length(bbbq::get_mhc1_haplotypes())

tibbles <- list()

for (i in seq_len(n_haplotypes)) {
  haplotype <- bbbq::get_mhc1_haplotypes()[i]
  t <- tibble::tibble(
    epitope = create_n_mers(protein_sequence, n),
    topology = create_n_mers(topology, n)
  )
  t$haplotype <- haplotype
  ic50_threshold <- epiprepreds::get_ic50_threshold(
    peptide_length = n,
    haplotype_name = haplotype,
    percentile = 0.02
  )
  t$ic50 <- bbbq::predict_ic50s(
    protein_sequence = protein_sequence,
    peptide_length = n,
    haplotype = haplotype,
    ic50_prediction_tool = "EpitopePrediction"
  )$ic50
  t$is_binder <- t$ic50 < ic50_threshold
  tibbles[[i]] <- t
}

t <- dplyr::bind_rows(tibbles)

t$is_tmh <- stringr::str_detect(t$topology, "1")
t$is_tmh_epitope <- t$is_binder & t$is_tmh

f_tmh <- mean(t$is_tmh)

t_perc <- t %>% dplyr::group_by(haplotype) %>% dplyr::summarise(f = mean(is_tmh_epitope) / mean(is_binder), f_tmh = mean(is_tmh), .groups = "keep")
library(ggplot2)
ggplot(t_perc, aes(x = haplotype, y = f)) +
  geom_col(col = "grey", fill = "grey") +
  geom_hline(yintercept = f_tmh, col = "red") +
  scale_y_continuous(
    "% epitopes overlapping\nwith transmembrane helix",
    labels = scales::percent
  )+
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.background = element_rect(fill = "white", color = "white")
  ) + ggsave("~/smimm11_perc_epitopes_tmh.png", width = 7, height = 7)
t_perc
