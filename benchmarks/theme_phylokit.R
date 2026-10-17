# House ggplot2 theme + palette for phylokit-mcp figures.
#
# No .R file existed anywhere in this repo before this one (checked with
# `git ls-files '*.R'` on 2026-10-16), so it is created here rather than
# reused. Modelled on plant-genomics-mcp's examples/arf_family/theme_pgmcp.R so
# the sibling servers' benchmark figures read alike.
#
# All ggplot2 styling lives here, sourced by the plotting scripts; a script
# adds only the geoms and scales specific to its figure. Restyle here, never
# inline.

theme_phylokit <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0, "pt"),
      legend.box.spacing = grid::unit(6, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title.position = "plot"
    )
}

# Agent benchmark verdicts (benchmarks/agent_qa/score.py). Right and wrong are
# green and red, declined the neutral grey; no_answer (no ANSWER line at all)
# gets the orange so a format failure is never read as either neighbour.
pal_verdict <- c(
  right = "#2E7D5B",
  wrong = "#C0392B",
  declined = "#BDBDBD",
  no_answer = "#E69F00"
)

# Per-run cost/time points: one neutral colour for every setup. Setups are
# named on the axis; a colour per setup would reuse green and grey, which mean
# right and declined in the verdict figures beside it.
phylokit_point_fill <- "grey35"
phylokit_point_outline <- "white"
phylokit_point_size <- 1.8
phylokit_point_stroke <- 0.5
