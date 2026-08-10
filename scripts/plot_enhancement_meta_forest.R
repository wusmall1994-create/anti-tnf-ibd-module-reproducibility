options(stringsAsFactors = FALSE)

library(ggplot2)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_dir <- file.path(root, "results", "enhancement_gate_20260805")
out_dir <- file.path(data_dir, "figure5_meta")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

effects <- read.csv(file.path(data_dir, "cross_cohort_signature_effects.csv"), check.names = FALSE)
meta <- read.csv(file.path(data_dir, "random_effects_meta_analysis.csv"), check.names = FALSE)

modules <- c(
  "myeloid_inflammation",
  "epithelial_barrier_maturity",
  "epithelial_ifn_damage",
  "stromal_fibroinflammatory"
)
module_labels <- c(
  myeloid_inflammation = "a  Myeloid inflammation",
  epithelial_barrier_maturity = "b  Epithelial barrier maturity",
  epithelial_ifn_damage = "c  Epithelial IFN/damage",
  stromal_fibroinflammatory = "d  Stromal fibroinflammatory"
)
stratum_labels <- c(
  GSE16879_CD_colon = "GSE16879 CD colon",
  GSE16879_UC_colon = "GSE16879 UC colon",
  GSE23597_UC_colon = "GSE23597 UC colon",
  GSE92415_UC_colon = "GSE92415 UC colon",
  `E-MTAB-7604_IBD_colon` = "E-MTAB-7604 colon",
  `E-MTAB-7604_IBD_ileum` = "E-MTAB-7604 ileum"
)

d <- effects[effects$signature %in% modules, ]
d$label <- unname(stratum_labels[d$stratum])
d$label <- paste0(d$label, " (", d$n_responder, "/", d$n_non_responder, ")")
d$estimate <- d$hedges_g_R_minus_NR
d$low <- d$ci95_low_g
d$high <- d$ci95_high_g
d$type <- "Cohort stratum"

m <- meta[meta$signature %in% modules & meta$meta_scope == "four_independent_cohorts_primary", ]
m$label <- "REML pooled (4 cohorts)"
m$estimate <- m$pooled_g_reml
m$low <- m$ci95_low_hk
m$high <- m$ci95_high_hk
m$type <- "Pooled"
m$n_responder <- NA
m$n_non_responder <- NA

plot_data <- rbind(
  d[, c("signature", "label", "estimate", "low", "high", "type")],
  m[, c("signature", "label", "estimate", "low", "high", "type")]
)

label_order <- c(
  "GSE16879 CD colon (12/7)",
  "GSE16879 UC colon (8/16)",
  "GSE23597 UC colon (25/7)",
  "GSE92415 UC colon (32/27)",
  "E-MTAB-7604 colon (10/17)",
  "E-MTAB-7604 ileum (9/8)",
  "REML pooled (4 cohorts)"
)
plot_data$label <- factor(plot_data$label, levels = rev(label_order))
plot_data$signature <- factor(plot_data$signature, levels = modules, labels = module_labels[modules])

hetero <- m[, c("signature", "I2_percent", "Q_p")]
hetero$signature <- factor(hetero$signature, levels = modules, labels = module_labels[modules])
hetero$annotation <- sprintf("I² = %.1f%%; Q p = %.3f", hetero$I2_percent, hetero$Q_p)
hetero$x <- -3.05
hetero$y <- 0.58

p <- ggplot(plot_data, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777", linetype = "dashed") +
  geom_errorbarh(aes(xmin = low, xmax = high, colour = type), height = 0, linewidth = 0.55) +
  geom_point(aes(fill = type, shape = type, size = type), colour = "white", stroke = 0.45) +
  facet_wrap(~signature, ncol = 2) +
  geom_text(data = hetero, aes(x = x, y = y, label = annotation), inherit.aes = FALSE,
            hjust = 0, vjust = 0, size = 2.35, colour = "#4C4C4C") +
  scale_colour_manual(values = c("Cohort stratum" = "#6E9BC4", "Pooled" = "#173F5F")) +
  scale_fill_manual(values = c("Cohort stratum" = "#8FB9DD", "Pooled" = "#173F5F")) +
  scale_shape_manual(values = c("Cohort stratum" = 21, "Pooled" = 23)) +
  scale_size_manual(values = c("Cohort stratum" = 2.5, "Pooled" = 3.6)) +
  scale_x_continuous(limits = c(-3.1, 3.1), breaks = seq(-3, 3, 1)) +
  labs(x = "Hedges g (responder minus non-responder)", y = NULL) +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.4, colour = "#222222"),
    axis.text.x = element_text(size = 6.4, colour = "#222222"),
    axis.title.x = element_text(size = 7.2, margin = margin(t = 5)),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.6, face = "bold", hjust = 0),
    panel.spacing.x = unit(10, "pt"),
    panel.spacing.y = unit(12, "pt"),
    legend.position = "none",
    plot.margin = margin(6, 8, 6, 6)
  )

write.csv(plot_data, file.path(out_dir, "Figure_5_source_data.csv"), row.names = FALSE)

width_in <- 183 / 25.4
height_in <- 145 / 25.4
grDevices::tiff(file.path(out_dir, "Figure_5.tiff"), width = width_in, height = height_in,
                units = "in", res = 600, compression = "lzw", type = "cairo", family = "Arial")
print(p)
dev.off()

grDevices::cairo_pdf(file.path(out_dir, "Figure_5.pdf"), width = width_in, height = height_in,
                     family = "Arial")
print(p)
dev.off()

grDevices::svg(file.path(out_dir, "Figure_5.svg"), width = width_in, height = height_in,
               family = "Arial")
print(p)
dev.off()

grDevices::png(file.path(out_dir, "Figure_5_preview.png"), width = width_in, height = height_in,
               units = "in", res = 200, type = "cairo")
print(p)
dev.off()
