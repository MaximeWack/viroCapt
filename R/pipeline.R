#' Plot depth
#'
#' @export
#' @param stem Stem name for sam file and plot file
#' @param limit Limit the number of genotypes to display
plot_depth <- function(stem, limit = 5)
{
  stringr::str_c(stem, ".rds") %>%
    readRDS -> sam

    sam %>%
      limit_genotypes(limit) %>%
      ggplot_depth -> Plot

    stringr::str_c(stem, ".png") %>%
      ggplot2::ggsave(Plot)
}


#' Create the profile object
#'
#' @export
#' @param sam Stem name of the sam file
create_profile <- function(stem)
{
  stringr::str_c(stem, ".sam") %>%
    read_sam -> sam

  if (length(sam) > 0)
  {
    sam %>%
      dplyr::mutate(parsed = cigar %>% parse_cigar) %>%
      tidyr::unnest() %>%
      read_depth %>%
      downsample %>%
      saveRDS(stringr::str_c(stem, ".rds"))
  } else
  {
    tibble::tibble(genotype = factor("No result"), pos = 0, n = 0) %>%
      saveRDS(stringr::str_c(stem, ".rds"))
  }
}

#' Create the visualisation object
#'
#' @export
#' @param sam Stem name of the rds file for the sequencing depth
#' @param summ Stem name of the blat summary file
create_visu <- function(sam, summ, visu)
{
  stringr::str_c(summ, "_summary.tsv") %>%
    read_summary -> summ_blat

  stringr::str_c(sam, ".rds") %>%
    readRDS -> sam

  list(summary = summ_blat,
       profile = sam) %>%
  saveRDS(visu)
}


#' Plot the final plot
#'
#' @export
#' @param sam Stem name of the sam file
#' @param summ Stem name of the blat summary file
#' @param limit Limit the number of genotypes to display
plot_final <- function(sam, summ, limit = 1)
{
  stringr::str_c(summ, "_summary.tsv") %>%
    read_summary -> summ_blat

  stringr::str_c(sam, ".rds") %>%
    readRDS -> sam

  if (length(summ_blat) > 0)
  {
    sam %>%
      limit_genotypes(limit) -> sam

    summ_blat %>%
      dplyr::semi_join(sam, by = "genotype") -> summ_blat

    sam %>%
      ggplot_final(summ_blat) -> Plot

    stringr::str_c(summ, ".png") %>%
      ggplot2::ggsave(Plot)
  } else
  {
    cat(NULL, file = stringr::str_c(summ, ".png"))
  }
}


#' Extract non-viral sequences as fasta files
#'
#' @export
#' @param samfile A local alignment sam file
#' @param fastafile The fasta file to write the result to
extract_fasta <- function(samfile, fastafile)
{
  samfile %>%
    read_sam -> sam

  if (length(sam) > 0)
  {
    sam %>%
      dplyr::filter(cigar %>% stringr::str_detect("S")) %>%
      dplyr::mutate(parsed = cigar %>% parse_cigar,
                    parsed = purrr::map2(parsed, pos, extract_features)) %>%
      tidyr::unnest() %>%
      dplyr::filter(length_read > 25,
                    feature %in% c("left", "right")) %>%
      extract_unaligned %>%
      write_fasta(fastafile)
  } else
  {
    cat(NULL, file = fastafile)
  }
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
    read_blat -> blat

  if (length(blat) > 0)
  {
    blat %>%
      clean_blat %>%
      tag_blat(fasta) %>%
      write_blat(cleaned)
  } else
  {
    cat(NULL, file = cleaned)
  }
}


#' Summarise blat file
#'
#' @export
#' @param blatfile A headerless cleaned blat file name
#' @param summarised A file name for output
summarised_blat <- function(blatfile, summarised)
{
  blatfile %>%
    readr::read_tsv() -> blat

  if (length(blat) > 0)
  {
    blat %>%
      summarise_blat %>%
      write_blat(summarised)
  } else
  {
    cat(NULL, file = summarised)
  }
}
