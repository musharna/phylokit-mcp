# Figures for the agent benchmark, from results.jsonl and from summary.tsv,
# which analyze.py writes (run it first):
#
#   uv run python benchmarks/agent_qa/analyze.py
#   Rscript benchmarks/agent_qa/plot.R [results.jsonl[,results_read.jsonl]] [outdir]
#
# Four figures:
#   agent_qa.png             right / wrong / declined / invalid per setup (the README one)
#   agent_qa_accuracy.png    accuracy with 95% intervals: setups and no-LLM baselines
#   agent_qa_cost_time.png   cost and seconds per answer, every run
#   agent_qa_by_category.png verdicts per setup x kind of question
# Verdict bars count answers, not questions: questions x runs per setup.
# Every number on the figures is computed here from the log; nothing is typed
# in by hand. Green, red and grey mean right, wrong and declined in every
# figure; setups are told apart by their axis labels, not by colour.

library(ggplot2)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
results_path <- if (length(args) >= 1) args[1] else "benchmarks/agent_qa/results.jsonl"
outdir <- if (length(args) >= 2) args[2] else "benchmarks/agent_qa"
source("benchmarks/theme_phylokit.R")

# One results file or several, comma-separated, as analyze.py takes them.
r <- do.call(rbind, lapply(strsplit(results_path, ",")[[1]], function(p) {
  x <- jsonlite::stream_in(file(p), verbose = FALSE)
  x[, c("arm", "qid", "category", "rep", "verdict", "valid", "wrong_id", "cost_usd", "seconds", "tool_calls")]
}))
stopifnot(!anyDuplicated(r[, c("arm", "qid", "rep")]))
stopifnot(all(r$verdict %in% names(pal_verdict)), "valid" %in% names(r))
# One rule for every setup: an invalid row is drawn as invalid, whatever its
# reply scored.
r$verdict[!r$valid] <- "invalid"

all_arm_labels <- c(
  bare = "Claude alone",
  read = "Claude + Read only (control)",
  alt = "Claude + Python, MAFFT, web",
  mcp = "Claude + phylokit-mcp"
)
stopifnot(all(r$arm %in% names(all_arm_labels)))
arm_labels <- all_arm_labels[names(all_arm_labels) %in% r$arm]
r$arm_label <- factor(arm_labels[r$arm], levels = rev(arm_labels))
r$verdict <- factor(r$verdict, levels = names(pal_verdict))

cat_labels <- c(
  sister_taxon = "Sister taxon in the ML tree",
  rf_distance = "Robinson-Foulds distance",
  homologous_position = "Homologous residue after alignment",
  model_selection = "Model with the lowest BIC",
  pairwise_distance = "JC69 distance between two taxa",
  informative_sites = "Parsimony-informative sites"
)
stopifnot(all(r$category %in% names(cat_labels)))

n_per_arm <- table(r$arm)
stopifnot(length(unique(n_per_arm)) == 1)
n_ans <- as.numeric(n_per_arm[1])
n_q <- length(unique(r$qid))
n_rep <- length(unique(r$rep))

# The tables come from analyze.py; this script only draws. It checks that
# summary.tsv was written from these same results.
summary_path <- file.path(outdir, "summary.tsv")
if (!file.exists(summary_path)) stop("no ", summary_path, ": run analyze.py first")
s <- read.delim(summary_path, na.strings = "NA")
sa <- s[s$level == "arm", ]
stopifnot(setequal(sa$arm, names(arm_labels)))
for (a in sa$arm) {
  stopifnot(sa$n[sa$arm == a] == sum(r$arm == a), sa$right[sa$arm == a] == sum(r$arm == a & r$verdict == "right"))
}

setup_note <- sprintf("Claude Opus 5.5 in Claude Code, %d questions × %d runs.", n_q, n_rep)
key_note <- paste0(
  "Answers fixed before any run by the generator: trees drawn at random, sequences evolved\n",
  "down them; the key never runs the server's code."
)
# Count labels: white on green and red, dark on the light grey of declined.
verdict_text <- c(`TRUE` = "grey15", `FALSE` = "white")
legend_counts <- function(v) {
  sprintf("%s (%s)", sub("_", " ", v), vapply(v, function(k) {
    paste(vapply(names(arm_labels), function(a) {
      as.character(sum(r$verdict == k & r$arm == a))
    }, ""), collapse = " / ")
  }, ""))
}

# Figure 1: one bar per setup, segments = verdict share, labelled with counts.
counts <- as.data.frame(table(arm_label = r$arm_label, verdict = r$verdict))
counts <- counts[counts$Freq > 0, ]
counts$share <- counts$Freq / n_ans
p1 <- ggplot(counts, aes(x = share, y = arm_label, fill = verdict)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = ifelse(share >= 0.05, Freq, ""), colour = verdict == "declined"),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    fontface = "bold", size = 4.2, show.legend = FALSE
  ) +
  scale_colour_manual(values = verdict_text, guide = "none") +
  # Slices too thin to hold a number are counted in the legend instead, per
  # setup in the order the title names.
  scale_fill_manual(
    values = pal_verdict, labels = legend_counts, drop = TRUE,
    name = paste0("Answers (", paste(sub("Claude \\+ |Claude ", "", arm_labels), collapse = " / "), "):")
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = sprintf("Share of %d answers per setup", n_ans), y = NULL,
    title = "Phylogenetics questions: right, wrong, declined or invalid",
    subtitle = paste0(setup_note, "\n", key_note)
  ) +
  theme_phylokit(base_size = 13) +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(8, 24, 8, 8))
ggsave(file.path(outdir, "agent_qa.png"), p1, width = 11, height = 1.8 + 0.7 * length(arm_labels), dpi = 200, bg = "white")

# Figure 2: accuracy with bootstrap 95% intervals over questions, Claude setups
# (filled) beside the no-LLM baselines (hollow), each on the questions it
# applies to.
# The subset rows (analyze.py NEEDS_COMPUTE) are the kinds the Read-only
# control did not solve by hand, chosen after it ran: labelled post hoc.
acc <- s[s$level %in% c("arm", "arm_x_subset", "baseline"), ]
acc <- acc[order(match(acc$level, c("arm", "arm_x_subset", "baseline"))), ]
acc$label <- ifelse(
  acc$level == "arm", sprintf("%s (all %d q.)", unname(all_arm_labels[acc$arm]), acc$questions),
  ifelse(acc$level == "arm_x_subset",
    sprintf("%s (sister, model, informative sites: %d q., post hoc)", unname(all_arm_labels[acc$arm]), acc$questions),
    sprintf("Baseline: %s (%s, %d q.)", acc$arm,
      ifelse(acc$category == "all", "all kinds", unname(cat_labels[acc$category])), acc$questions)
  )
)
stopifnot(!anyDuplicated(acc$label))
acc$label <- factor(acc$label, levels = rev(acc$label))
acc$kind <- ifelse(acc$level == "baseline", "No-LLM baseline (one answer per question)", "Claude setup (all answers)")
p_acc <- ggplot(acc, aes(x = accuracy, y = label)) +
  geom_errorbar(aes(xmin = accuracy_ci_lo, xmax = accuracy_ci_hi), orientation = "y", width = 0.25, colour = phylokit_ink) +
  geom_point(aes(shape = kind), size = 3, fill = phylokit_ink, colour = phylokit_ink, stroke = 0.9) +
  scale_shape_manual(values = c(21, 1), name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Share right (invalid rows count as not right).\nBars: bootstrap 95% interval over questions",
    y = NULL, title = "Accuracy of each setup beside what a fixed guess or a shortcut gets",
    subtitle = paste0(
      setup_note, "\nBaselines are scored by the same scorer on the same data files; the best constant ",
      "is chosen with the key, an upper bound."
    )
  ) +
  theme_phylokit(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(8, 24, 8, 8))
ggsave(file.path(outdir, "agent_qa_accuracy.png"), p_acc, width = 13.5, height = 1.8 + 0.42 * nrow(acc), dpi = 200, bg = "white")

# Figure 3: cost and seconds per answer, one point per run, box over them.
long <- rbind(
  data.frame(arm_label = r$arm_label, metric = "Cost per answer (US$, log scale)", value = r$cost_usd),
  data.frame(arm_label = r$arm_label, metric = "Seconds per answer (log scale)", value = r$seconds)
)
p2 <- ggplot(long, aes(x = value, y = arm_label)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "white", colour = "grey20") +
  geom_point(
    position = position_jitter(height = 0.18, seed = 1),
    shape = 21, size = phylokit_point_size, stroke = phylokit_point_stroke,
    fill = phylokit_point_fill, colour = phylokit_point_outline, alpha = 0.55
  ) +
  facet_wrap(~metric, ncol = 2, scales = "free_x") +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  labs(
    x = NULL, y = NULL,
    title = "Cost and time per answer",
    subtitle = paste0(
      setup_note, "\nOne point per run; the box spans the middle half, the line in it is the median."
    )
  ) +
  theme_phylokit(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
    plot.margin = margin(8, 24, 8, 8)
  )
ggsave(file.path(outdir, "agent_qa_cost_time.png"), p2, width = 12, height = 1.8 + 0.7 * length(arm_labels), dpi = 200, bg = "white")

# Figure 4: verdicts per setup x kind of question; strips name the counts.
n_cat <- tapply(r$qid, r$category, function(x) length(unique(x)))
strip <- setNames(
  sprintf("%s (%d questions × %d runs)", cat_labels, as.integer(n_cat[names(cat_labels)]), n_rep),
  names(cat_labels)
)
r$cat_strip <- factor(strip[r$category], levels = strip)
cc <- as.data.frame(table(arm_label = r$arm_label, cat_strip = r$cat_strip, verdict = r$verdict))
cc <- cc[cc$Freq > 0, ]
p3 <- ggplot(cc, aes(x = Freq, y = arm_label, fill = verdict)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = Freq, colour = verdict == "declined"),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    fontface = "bold", size = 3.4, show.legend = FALSE
  ) +
  scale_colour_manual(values = verdict_text, guide = "none") +
  facet_wrap(~cat_strip, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = pal_verdict, name = NULL, drop = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    x = "Answers", y = NULL, title = "Right, wrong, declined or invalid, by kind of question",
    subtitle = setup_note
  ) +
  theme_phylokit(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
    plot.margin = margin(8, 24, 8, 8)
  )
ggsave(file.path(outdir, "agent_qa_by_category.png"), p3, width = 14, height = 8, dpi = 200, bg = "white")

print(sa[, c("arm", "n", "right", "wrong", "declined", "invalid", "accuracy", "accuracy_ci_lo", "accuracy_ci_hi")])
