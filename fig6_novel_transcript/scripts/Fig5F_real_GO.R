# Panel F: GO Term Enrichment from REAL novel transcript annotation
# Parses GO column from new_transcripts.all_annotation.xls (14,381 transcripts)
library(data.table); library(ggplot2)
source("../../shared/palette/color_palette.R")
source("fig5_muted_palette.R")
out_dir <- "../output"

# ---- Parse GO terms from annotation XLS ----
annot_file <- "../data/new_transcripts.all_annotation.xls"
if (!file.exists(annot_file)) stop("Annotation file not found: ", annot_file)

cat("Parsing GO terms from annotation file...\n")
annot <- fread(annot_file, sep="\t", header=TRUE, quote="", showProgress=FALSE)
cat(sprintf("  Total transcripts: %d\n", nrow(annot)))

go_col <- annot$GO[annot$GO != "" & annot$GO != "None" & !is.na(annot$GO)]
cat(sprintf("  With GO annotations: %d\n", length(go_col)))

# Parse GO entries: "GO:XXXXXXX//TERM_NAME`category;..."
all_go <- unlist(strsplit(go_col, ";"))
all_go <- trimws(all_go[all_go != ""])

parse_one_go <- function(x) {
  parts <- strsplit(x, "//")[[1]]
  if (length(parts) < 2) return(NULL)
  term_parts <- strsplit(parts[2], "`")[[1]]
  if (length(term_parts) < 2) return(NULL)
  list(term = term_parts[1], category = term_parts[2])
}

parsed <- lapply(all_go, parse_one_go)
parsed <- parsed[!sapply(parsed, is.null)]
go_dt <- rbindlist(lapply(parsed, as.data.table))
cat(sprintf("  Parsed GO entries: %d, unique terms: %d\n", nrow(go_dt), uniqueN(go_dt$term)))

# Count and rank
go_counts <- go_dt[, .N, by = .(term, category)][order(-N)]

# Map to simplified categories
go_counts[, Category := fcase(
  grepl("cellular_component", category), "Cellular Component",
  grepl("molecular_function", category), "Molecular Function",
  grepl("biological_process", category), "Biological Process",
  default = "Other"
)]

# Top 30 non-Other terms
top_go <- go_counts[Category != "Other"][1:30]
top_go[, Pct := N / nrow(annot) * 100]
top_go[, GO_Term := factor(term, levels = rev(term))]

# ---- Abbreviate long GO terms (2026-06-23) ----
go_abbrev <- c(
  "regulation of transcription by RNA polymerase II" = "RNA Pol II transcription regulation",
  "positive regulation of transcription by RNA polymerase II" = "RNA Pol II transcription activation",
  "negative regulation of transcription by RNA polymerase II" = "RNA Pol II transcription repression",
  "RNA polymerase II cis-regulatory region sequence-specific DNA binding" = "RNA Pol II cis-regulatory DNA binding",
  "DNA-binding transcription factor activity, RNA polymerase II-specific" = "RNA Pol II-specific TF activity",
  "endoplasmic reticulum membrane" = "ER membrane",
  "endoplasmic reticulum" = "ER",
  "perinuclear region of cytoplasm" = "perinuclear cytoplasm",
  "mitochondrial inner membrane" = "mitochondrial inner membrane",
  "extracellular space" = "extracellular space",
  "identical protein binding" = "identical protein binding"
)
top_go[, term_short := ifelse(term %in% names(go_abbrev), go_abbrev[term], term)]
top_go[, GO_Term := factor(term_short, levels = rev(unique(term_short[order(match(term, rev(levels(GO_Term))))])))]

cat(sprintf("  Top 30 GO terms selected for plotting\n"))

# ---- Plot ----
go_cat_colors <- c(
  "Cellular Component" = "#8CBFD0",
  "Molecular Function" = "#D55E5E",
  "Biological Process" = "#9C78B8"
)

p_f <- ggplot(top_go, aes(x = N, y = GO_Term, fill = Category)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", N, Pct)), hjust = -0.1, size = 3.2) +
  scale_fill_manual(values = go_cat_colors, name = "GO Category") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, max(top_go$N) * 1.25)) +
  labs(x = "Number of Novel Transcript Host Genes", y = NULL,
       title = "Top GO Terms of Novel Transcript Host Genes",
       subtitle = sprintf("Parsed from %d novel transcripts (%.0f%% have GO annotations); %d unique GO terms",
                          nrow(annot), length(go_col)/nrow(annot)*100, uniqueN(go_dt$term))) +
  theme_drs + theme(
    panel.grid.major.y = element_blank(),
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill = "white", color = "grey80")
  )

pdf(file.path(out_dir, "Fig6D_GO_Enrichment.pdf"), w = 11, h = 8)
print(p_f)
dev.off()
png(file.path(out_dir, "Fig6D_GO_Enrichment.png"), w = 11, h = 8, units = "in", res = 300)
print(p_f)
dev.off()

# Save parsed GO data
fwrite(top_go, file.path(out_dir, "Fig6F_GO_data.csv"))
cat(sprintf("\nFig5F saved. %d GO terms parsed from real annotation.\n", nrow(top_go)))
