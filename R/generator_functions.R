
#' Generate a final plot function
#'
#' @export
#' @param write Write plot on disk?
#' @param fastqdir Directory of fastq files
#' @param dir Directory in which to write the plot if write == T
#' @return A function accepting a local plot object, a blat summary, and a fastq name, and returning a plot object
Gplot_final <- function(write = F, fastqdir = NULL, dir = NULL)
{
  if (write)
    dir.create(paste0(fastqdir, "/", dir), showWarnings = F, recursive = T)

  function(local_plot, summ_blat, fastq)
  {
    summ_blat %>% filter(!is.na(quality), quality >= "HC") -> summ_blat

    local_plot$data$n %>% max -> max_n
    local_plot +
      ggplot2::geom_vline(data = summ_blat, ggplot2::aes(xintercept = position, color = chr, alpha = quality)) -> p

    if (summ_blat %>% filter(feature == "left") %>% nrow > 0)
      p + ggplot2::geom_segment(data = summ_blat %>% filter(feature == "left"), ggplot2::aes(x = position, xend = position - 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    if (summ_blat %>% filter(feature == "right") %>% nrow > 0)
      p + ggplot2::geom_segment(data = summ_blat %>% filter(feature == "right"), ggplot2::aes(x = position, xend = position + 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    p +
      ggplot2::scale_alpha_ordinal(limits = c("HC", "HCT"), range = c(0.5, 1)) +
      ggplot2::theme_classic() -> Plot

    if (write)
      ggplot2::ggsave(plot = Plot, filename = paste0(fastqdir, "/", dir, "/", fastq, ".png"))
  }
}
