#' Plot depth
#'
#' @export
#' @param stem Stem name for sam file and plot file
#' @param Save Whether to save the data of the plot as an rds file
plot_depth <- function(stem, Save = F)
{
  stringr::str_c(stem, ".sam") %>%
    read_sam %>%
    dplyr::mutate(parsed = cigar %>% parse_cigar) %>%
    tidyr::unnest() %>%
    read_depth %>%
    # normalise_depth(qc_norm) %>%
    downsample %>%
    # MA(15) %>%
    ggplot_depth -> Plot

  stringr::str_c(stem, ".png") %>%
    ggplot2::ggsave(Plot)

    if (Save)
      saveRDS(Plot, stringr::str_c(stem, ".rds"))
}


#' Extract non-HPV sequences as fasta files
#'
#' @export
#' @param samfile A local alignment sam file
#' @param fastafile The fasta file to write the result to
extract_fasta <- function(samfile, fastafile)
{
  samfile %>%
    read_sam %>%
    dplyr::filter(cigar %>% stringr::str_detect("S")) %>%
    dplyr::mutate(parsed = cigar %>% parse_cigar,
                  parsed = purrr::map2(parsed, pos, extract_features)) %>%
    tidyr::unnest() %>%
    dplyr::filter(length_read > 25,
           feature %in% c("left", "right")) %>%
    extract_unaligned %>%
    write_fasta(fastafile)
}


#' Clean blat file
#'
#' @export
#' @param blatfile A headerless blat file name
#' @param fastafile A headerless fasta file name
#' @param cleaned A file name for output
cleaned_blat <- function(blatfile, fastafile, cleaned)
{
  fastafile %>%
    read_fasta -> fasta

  blatfile %>%
    read_blat %>%
    clean_blat %>%
    tag_blat(fasta) %>%
    write_blat(cleaned)
}


#' Summarise blat file
#'
#' @export
#' @param blatfile A headerless cleaned blat file name
#' @param summarised A file name for output
summarised_blat <- function(blatfile, summarised)
{
  blatfile %>%
    readr::read_tsv() %>%
    summarise_blat %>%
    write_blat(summarised)
}
