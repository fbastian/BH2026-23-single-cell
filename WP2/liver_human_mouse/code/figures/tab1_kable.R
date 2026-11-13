# Table 1 - SAMap matches mouse (Guilliams 2022 rebuilt) -> human (HLiCA v1.0), score >= 0.4, with top gene pairs.
# kableExtra 1.4.0 (env masld_tables): LaTeX (booktabs, longtable, grouped by lineage) -> PDF via tectonic (env latex); HTML copy.
# Input/Output: /mnt/masld_vault/cross_species/figures/tables/table1_samap_matches.{tsv,tex,pdf,html}
suppressPackageStartupMessages({library(kableExtra); library(dplyr)})
D <- "/mnt/masld_vault/cross_species/figures/tables"
t <- read.delim(file.path(D, "table1_samap_matches.tsv"), check.names = FALSE) %>% select(-`Gene pairs`)
lin <- t$Lineage; t <- t %>% select(-Lineage)
idx <- table(factor(lin, levels = unique(lin)))
cap <- paste0("SAMap cell-type matches between the rebuilt Guilliams et al. 2022 mouse liver atlas and HLiCA v1.0 (mapping score $\\geq$ 0.4; ", nrow(t), " matches). Top gene pairs: GenePairFinder (align\\_thr 0.1), human/mouse symbol; one name when both agree.")
tex <- kbl(t, format = "latex", booktabs = TRUE, longtable = TRUE, caption = cap, escape = TRUE, linesep = "") %>%
  kable_styling(latex_options = c("repeat_header", "striped"), font_size = 7, stripe_color = "gray!8") %>%
  column_spec(1, width = "2.3cm") %>% column_spec(2, width = "3.8cm") %>% column_spec(3, width = "1.4cm") %>%
  column_spec(4, width = "13.5cm", italic = TRUE) %>%
  pack_rows(index = idx, latex_gap_space = "0.3em") %>%
  collapse_rows(columns = 1, latex_hline = "none", valign = "top")
doc <- c("\\documentclass{article}", "\\usepackage[a4paper,landscape,margin=1.2cm]{geometry}",
         "\\usepackage[table]{xcolor}", "\\usepackage{booktabs,longtable,array,multirow,colortbl,makecell}",
         "\\usepackage{fontspec}", "\\setmainfont{Arial}", "\\setsansfont{Arial}", "\\pagestyle{empty}",
         "\\begin{document}", as.character(tex), "\\end{document}")
writeLines(doc, file.path(D, "table1_samap_matches.tex"))
h <- kbl(t, caption = gsub("\\$\\\\geq\\$", "&ge;", gsub("\\\\_", "_", cap)), escape = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Arial", font_size = 12) %>% column_spec(4, italic = TRUE) %>%
  pack_rows(index = idx) %>% collapse_rows(columns = 1, valign = "top")
save_kable(h, file.path(D, "table1_samap_matches.html"), self_contained = TRUE)
cat("written\n")
