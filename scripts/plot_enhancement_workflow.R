library(grid)

out_dir <- file.path("results", "enhancement_gate_20260805", "figure1_workflow")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

navy <- "#0B2A4A"
blue <- "#2D74B8"
teal <- "#2F8F8B"
orange <- "#C86F1A"
grey <- "#596273"
light_blue <- "#EAF3FB"
light_teal <- "#E8F5F2"
light_orange <- "#FFF2E5"
light_grey <- "#F4F6F8"

box <- function(x, y, w, h, fill, border, title, lines, title_col = navy) {
  grid.roundrect(x, y, w, h, r = unit(0.035, "snpc"),
                 gp = gpar(fill = fill, col = border, lwd = 2.5))
  grid.text(title, x = unit(x - w/2 + 0.025, "npc"), y = unit(y + h/2 - 0.045, "npc"),
            just = c("left", "top"), gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 16, col = title_col))
  grid.text(paste(lines, collapse = "\n"), x = unit(x - w/2 + 0.025, "npc"), y = unit(y + h/2 - 0.105, "npc"),
            just = c("left", "top"), gp = gpar(fontfamily = "sans", fontsize = 11.5, col = grey, lineheight = 1.25))
}

draw_arrow <- function(x0, y0, x1, y1, col = "#6B7280") {
  grid.lines(c(x0, x1), c(y0, y1), arrow = grid::arrow(type = "closed", length = unit(0.018, "npc")),
             gp = gpar(col = col, lwd = 2.5))
}

draw_workflow <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("Public-data workflow for cross-cohort anti-TNF mucosal modules",
            x = 0.045, y = 0.955, just = c("left", "top"),
            gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 25, col = navy))
  grid.text("Independent-cohort auditing, formal meta-analysis, covariate checks and single-cell contextualisation",
            x = 0.045, y = 0.905, just = c("left", "top"),
            gp = gpar(fontfamily = "sans", fontsize = 13, col = grey))

  box(0.18, 0.72, 0.27, 0.23, light_grey, "#AEB8C4", "1. Public resources",
      c("Bulk: GSE16879 and GSE23597", "GSE92415 and E-MTAB-7604", "Single cell: TAURUS / GSE282122", "Candidate series retained for audit"))
  box(0.50, 0.72, 0.30, 0.23, light_blue, blue, "2. Eligibility gate",
      c("Include: four independent pretreatment", "mucosal cohorts with outcome labels", "Exclude overlap: GSE14580, GSE12251", "Exclude uncertain reuse: GSE73661"))
  box(0.835, 0.72, 0.285, 0.23, light_teal, teal, "3. Parallel analyses",
      c("Effect sizes + REML meta-analysis", "Mayo/age; tissue/drug adjustment", "TAURUS lineage context", "Arijs and OSM benchmarks"))
  draw_arrow(0.315, 0.72, 0.345, 0.72)
  draw_arrow(0.65, 0.72, 0.692, 0.72)

  box(0.26, 0.42, 0.38, 0.23, light_teal, teal, "Bulk synthesis",
      c("Myeloid inflammation: lower in responders", "g = -0.99 (95% CI -1.67 to -0.32)", "Epithelial barrier maturity: higher", "g = 0.56 (95% CI 0.18 to 0.94)"))
  box(0.74, 0.42, 0.38, 0.23, light_blue, blue, "Cellular and benchmark context",
      c("Resident-like macrophage composition", "contextualises the bulk myeloid signal", "AUC comparable with published signatures", "No evidence of predictive superiority"))
  draw_arrow(0.84, 0.605, 0.74, 0.535, teal)
  draw_arrow(0.78, 0.605, 0.26, 0.535, teal)

  box(0.50, 0.155, 0.78, 0.21, light_orange, orange, "Bounded conclusion",
      c("Pretreatment mucosal myeloid and barrier states show reproducible anti-TNF outcome associations",
        "across public cohorts. Prospective, clinically harmonised validation is required before", "treatment-selection use."), title_col = navy)
  draw_arrow(0.26, 0.305, 0.42, 0.235, orange)
  draw_arrow(0.74, 0.305, 0.58, 0.235, orange)
}

tiff(file.path(out_dir, "Figure_1.tiff"), width = 3600, height = 2400, res = 300, compression = "lzw")
draw_workflow()
dev.off()

png(file.path(out_dir, "Figure_1_preview.png"), width = 1800, height = 1200, res = 150)
draw_workflow()
dev.off()

pdf(file.path(out_dir, "Figure_1.pdf"), width = 12, height = 8, useDingbats = FALSE)
draw_workflow()
dev.off()
