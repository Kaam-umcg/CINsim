# plotting values

# copy number colors
copy_num_cols <- c(
  "gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
  "lemonchiffon", "dodgerblue", "chartreuse4", "lightcoral", "aquamarine2"
)
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

# plotting theme
cinsim_theme <- function() {
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(colour = "black"),
      aspect.ratio = 1
    )
}
