#' Read a headerless sam file
#'
#' @param filename The SAM filename
#' @return A dataframe with columns (read, genotype, pos, cigar, seq)
read_sam <- function(filename)
{
  suppressWarnings({
  readr::read_tsv(filename,
           col_names = c("read", "flag", "genotype", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq"),
           col_types = "ciciicciic") -> sam

  if (length(sam) > 0)
    sam %>%
      dplyr::select(read, genotype, pos, cigar, seq) %>%
      dplyr::mutate(genotype = forcats::fct_infreq(genotype))
  })
}


#' Compute the sequencing depth
#'
#' @param sam sam file object (dataframe)
#' @return A sequencing depth object
read_depth <- function(sam)
{
  sam %>%
    dplyr::filter(type == "M") %>%
    dplyr::mutate(start = pos + start_hpv - 1,
                  end = pos + end_hpv - 1) %>%
    dplyr::mutate(base = purrr::map2(start, end, base::seq)) %>%
    dplyr::group_by(genotype) %>%
    dplyr::do(data = .$base %>%
              unlist %>%
              tibble::tibble(pos = .) %>%
              dplyr::count(pos) %>%
              tidyr::complete(pos = 1:max(pos, na.rm = T), fill = list(n = 0))) %>%
    tidyr::unnest() %>%
    dplyr::ungroup()
}


#' Normalise a sequencing depth object using a reference
#'
#' @param depth A sequencing depth object produced by read_depth
#' @param qc_norm A normalization map for the depth object
#' @return A normalized sequencing depth object
normalise_depth <- function(depth, qc_norm)
{
  depth %>%
    dplyr::left_join(qc_norm, by = c("genotype", "pos")) %>%
    dplyr::mutate(n = n.x / n.y) %>%
    dplyr::select(-n.x, -n.y)
}


#' Downsample and object of sequencing depth
#'
#' The sequencing depth object has a width corresponding to the length of the aligned genome (~8000bp for HPV)
#' This is way too large for plotting on a regular screen
#' This function allows the downsampling of a depth object to speed up other operations, such as plotting
#'
#' @param depth A sequencing depth object, produced by read_depth
#' @param ratio The downsampling ratio, 5 by default
#' @return A downsampled sequencing depth object
downsample <- function(depth, ratio = 5)
{
  depth %>%
    dplyr::mutate(pos = ( (pos / ratio) %>% ceiling) * ratio) %>%
    dplyr::group_by(genotype, pos) %>%
    dplyr::summarise(n = mean(n, na.rm = T)) %>%
    dplyr::ungroup()
}


#' Compute moving average
#'
#' @param depths A nucleotide depth object
#' @param window Half-size of the window
#' @return A smoothed depth object
MA <- function(depths, window)
{
  depths %>%
    dplyr::filter(n > 0) %>%
    split(.$genotype) %>%
    lapply(function(x)
           {
             x %>%
               dplyr::mutate(roll = pos %>% purrr::map_dbl(~mean(x$n[x$pos %in% seq(. - window, . + window)], na.rm = T))) %>%
               tidyr::complete(pos = 1:max(pos, na.rm = T), fill = list(genotype = x$genotype %>% unique, roll = 0, n = 0))
           }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(n = roll)
}


#' Limit the number of genotypes in a sample
#'
#' @param depths A nucleotide depth object
#' @param limit The number of genotypes to keep
#' @return A depth object with only the top n genotypes
limit_genotypes <- function(depths, limit)
{
  depths %>%
    dplyr::filter(genotype %in% base::levels(genotype)[1:limit]) %>%
    dplyr::mutate(genotype = genotype %>% forcats::fct_drop())
}


#' Plot the nucleotide depth
#'
#' @param depth A sequencing depth object
#' @return A nucleotide depth plot
ggplot_depth <- function(depth)
{
  depth %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = pos, y = n) +
    ggplot2::geom_area() +
    ggplot2::facet_grid(~genotype) +
    ggplot2::theme(strip.text = ggplot2::element_text(angle = -90))
}


#' Parse a cigar string into a dataframe with type, length, start and end position
#'
#' @param cigar vector of CIGAR string
#' @return parsed CIGAR as a dataframe
parse_cigar <- function(cigar)
{
  cigar %>%
    gregexpr(pattern = "\\d+(M|I|D|S|N|H|P|=|X)") %>%
    regmatches(x = cigar) %>%
    lapply(function(x)
               {
                 type <- substring(x, nchar(x))
                 length <- substring(x, 0, nchar(x) -1) %>% as.numeric

                 length_read <- length
                 length_read[type == "D"] <- 0

                 length_hpv <- length
                 length_hpv[type == "I"] <- 0

                 start_read <- cumsum(length_read) - length_read + 1
                 start_hpv <- cumsum(length_hpv) - length_hpv + 1

                 end_read <- cumsum(length_read)
                 end_hpv <- cumsum(length_hpv)

                 data.frame(length_read = length_read, length_hpv = length_hpv, type = type, start_read = start_read, end_read = end_read, start_hpv = start_hpv, end_hpv = end_hpv, stringsAsFactors = F)
               })
}


#' Extract relevant features from a parsed CIGAR string
#'
#' @param parsed_cigar A parsed CIGAR as a dataframe (as returned by parse_cigar)
#' @param pos The matching position on HPV for that CIGAR string
#' @return A dataframe of features with their position on HPV
extract_features <- function(parsed_cigar, pos)
{
  parsed_cigar$length_read %>% sum -> read_length

  parsed_cigar %>%
    dplyr::mutate(feature = dplyr::case_when(type == "S" & start_read == 1 ~ "left",
                                             type == "S" & end_read == read_length ~ "right",
                                             T ~ type),
                  feature_pos = dplyr::case_when(feature == "D" ~ pos + start_hpv - 1,
                                                 feature == "I" ~ pos + start_hpv- 2,
                                                 feature == "left" ~ pos + end_hpv,
                                                 feature == "right" ~ pos + start_hpv - 2)) -> df

  if (df$feature[1] == "left")
    df$feature_pos <- df$feature_pos - df$length_read[1]

  df
}


#' Extract unaligned sequences
#'
#' @param parsed A sam object (dataframe) with parsed CIGARs and extracted features
#' @return A sam object with unaligned sequences extracted
extract_unaligned <- function(parsed)
{
  parsed %>%
    dplyr::mutate(nalign_seq = substring(seq, start_read, end_read)) %>%
    dplyr::arrange(feature_pos, dplyr::desc(length_read)) %>%
    dplyr::select(genotype, read, feature, feature_pos, nalign_seq) %>%
    droplevels
}


#' Write a fasta file from a parsed sam file
#'
#' @export
#' @param parsed_sam A parsed sam object (cigar went through parse_cigar and extract_features)
#' @param fastafile A file to write to
write_fasta <- function(parsed_sam, fastafile)
{
  paste(paste0(">",
               parsed_sam$read, "|",
               parsed_sam$genotype, "|",
               parsed_sam$feature, "|",
               parsed_sam$feature_pos),
        parsed_sam$nalign_seq,
        sep = "\n",
        collapse = "\n") %>%
    cat(file = fastafile)
}


#' Read a fasta file
#'
#' @export
#' @param fastafile A headerless fasta file to read
#' @return A fasta file as a dataframe
read_fasta <- function(fastafile)
{
  fastafile %>%
    readr::read_lines() -> raw

    tibble::tibble(desc = stringr::str_replace(raw[c(T, F)], "^>", ""),
           nalign_seq = raw[c(F, T)]) %>%
      tidyr::separate(desc, into = c("read", "genotype", "feature", "feature_pos"), sep = "\\|") %>%
      dplyr::mutate(feature_pos = feature_pos %>% as.numeric)
}


#' Write a blat file
#'
#' Copy of readr::write_tsv
#'
#' @export
write_blat <- readr::write_tsv


#' Read a headerless blat file
#'
#' @export
#' @param blatfile The file to read
#' @return A blatfile organized as a dataframe
read_blat <- function(blatfile)
{
  blatfile %>%
    readr::read_tsv(col_names = c("match","mis-match","rep.match","N's","Q gap count","Q gap bases","T gap count","T gap bases","strand","Q name","Q size","Q start","Q end","T name","T size","T start","T end","block count","blockSizes","qStarts","tStarts")) -> blat

  if (length(blat) > 0)
    blat %>%
      tidyr::separate(`Q name`, c("read", "genotype", "feature", "position"), sep = "\\|")
}


#' Clean a raw blat file
#'
#' @param blat A blat object (dataframe obtained from Gread_blat)
#' @return A cleaned, deduplicated blat object
clean_blat <- function(blat)
{
  blat %>%
    dplyr::mutate(prop = (`Q end` - `Q start`) / `Q size`,
                  chr_position = dplyr::case_when(feature == "left" & strand == "-" ~ `T start`,
                                                  feature == "left" & strand == "+" ~ `T end`,
                                                  feature == "right" & strand == "+" ~ `T start`,
                                                  feature == "right" & strand == "-" ~ `T end`)) %>%
  dplyr::filter(`Q gap bases` <= 5,
                `T gap bases` <= 5,
                prop > .9) %>%
  dplyr::rename(chr = `T name`, size = `Q size`) %>%
  dplyr::group_by(read, strand, genotype, feature, position, chr, chr_position, size) %>%
  dplyr::summarise(match = max(match, na.rm = T)) %>%
  dplyr::filter(size == max(size, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
}


#' Remove homologs (reads with multiple matches in different places)
#'
#' @param blat A clean blat object (dataframe)
#' @return A filtered blat object with homologs removed
remove_homologs <- function(blat)
{
  blat %>%
    dplyr::group_by(read, genotype, feature) %>%
    dplyr::add_count() %>%
    dplyr::filter(n == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n)
}


#' Find positions with more than one matching read
#'
#' @param blat A clean blat object (dataframe)
#' @return A filtered blat object with only concordant positions kept
find_concordant_positions <- function(blat)
{
  blat %>%
    dplyr::group_by(genotype, feature, position, chr, chr_position) %>%
    dplyr::add_count() %>%
    dplyr::filter(n > 1) %>%
    dplyr::arrange(chr, position, chr_position) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup()
}


#' Find positions with at least one left and right clip on the same chromosome
#'
#' @param blat A clean blat object (dataframe)
#' @return A filtered blat object with only two-sided clips
filter_two_sides <- function(blat)
{
  blat %>%
    dplyr::semi_join(blat %>%
                     dplyr::group_by(genotype, chr) %>%
                     dplyr::summarise(feat = dplyr::n_distinct(feature)) %>%
                     dplyr::filter(feat == 2) %>%
                     dplyr::ungroup())
}


#' Tag a blat object with quality measures
#'
#' @param blat A blat object (dataframe)
#' @param fasta The corresponding fasta file that produced the blat object
#' @return An tagged blat file enriched with the source sequence
tag_blat <- function(blat, fasta)
{
  compose_combn(H = remove_homologs,
                C = find_concordant_positions,
                T = filter_two_sides) -> comb_funs

  lapply(comb_funs, function(x, y) x(y), blat) -> res

  names(res) %>%
    lapply(function(x) {dplyr::mutate(res[[x]], quality = x)}) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(quality = quality %>% ordered(levels = c("T", "C", "H", "HT", "CT", "HC", "HCT"))) %>%
    dplyr::group_by(read, feature) %>%
    dplyr::filter(quality == max(quality, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::full_join(blat) %>%
    dplyr::mutate(position = position %>% as.numeric) %>%
    dplyr::left_join(fasta, by = c("genotype", "read", "feature", "position" = "feature_pos")) %>%
    dplyr::arrange(dplyr::desc(quality), chr, position, dplyr::desc(match)) %>%
    dplyr::filter(nchar(nalign_seq) == size) %>%
    dplyr::select(read, strand, sequence = nalign_seq, genotype, feature, position, chr, chr_position, match, quality)
}


#' Summarise the results of a tagged blat file
#'
#' @param blat A tagged blat file
#' @return A summary table of results
summarise_blat <- function(blat)
{
  blat %>%
    dplyr::mutate(quality = quality %>% ordered(levels = c("T", "C", "H", "HT", "CT", "HC", "HCT"))) %>%
    dplyr::group_by(genotype, feature, position, chr, chr_position) %>%
    dplyr::add_count() %>%
    dplyr::summarise(n = max(n, na.rm = T), quality = max(quality, na.rm = T), match = max(match, na.rm = T)) %>%
    dplyr::arrange(dplyr::desc(quality), dplyr::desc(n), dplyr::desc(match)) %>%
    dplyr::ungroup()
}


#' Read a blat summary file
#'
#' @param summ_blat Name of the summary file
#' @return A blat summary object
read_summary <- function(summ_blat)
{
  readr::read_tsv(summ_blat) -> blat

  if (length(blat) > 0)
    blat %>%
      dplyr::mutate(quality = quality %>% ordered(levels = c("T", "C", "H", "HT", "CT", "HC", "HCT")))
}


#' Plot the final plot with insertions
#'
#' @param depth A sequencing depth object
#' @param summ_blat A blat summary object
#' @return A nucleotide depth plot
ggplot_final <- function(depth, summ_blat)
{
  depth %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = pos, y = n) +
    ggplot2::geom_area() +
    ggplot2::facet_grid(~genotype) +
    ggplot2::theme(strip.text = ggplot2::element_text(angle = -90)) -> local_plot

    summ_blat %>% dplyr::filter(!is.na(quality), quality >= "HC") -> summ_blat

    max(depth$n) -> max_n

    local_plot +
      ggplot2::geom_vline(data = summ_blat, ggplot2::aes(xintercept = position, color = chr, alpha = quality)) -> p

    if (summ_blat %>% dplyr::filter(feature == "left") %>% nrow > 0)
      p + ggplot2::geom_segment(data = summ_blat %>% dplyr::filter(feature == "left"), ggplot2::aes(x = position, xend = position - 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    if (summ_blat %>% dplyr::filter(feature == "right") %>% nrow > 0)
      p + ggplot2::geom_segment(data = summ_blat %>% dplyr::filter(feature == "right"), ggplot2::aes(x = position, xend = position + 100, y = -max_n/20, yend = -max_n/20, color = chr, alpha = quality)) -> p

    p +
      ggplot2::scale_alpha_ordinal(limits = c("HC", "HCT"), range = c(0.5, 1)) +
      ggplot2::theme_classic()
}
