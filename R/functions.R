#' Parse a raw sam output to a dataframe
#'
#' @param rawsam Raw sam file (character vector)
#' @return sam file object (dataframe)
parse_sam <- function(rawsam)
{
  data.frame(rawsam) %>%
    tidyr::separate(rawsam, sep = "\t", into = c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual")) %>%
    dplyr::mutate(flag = flag %>% as.numeric,
                  pos = pos %>% as.numeric,
                  mapq = mapq  %>% as.numeric,
                  pnext = pnext  %>% as.numeric,
                  tlen = tlen  %>% as.numeric,
                  rname = rname %>% factor)
}


#' Compute the nucleotide depth
#'
#' @export
#' @param sam sam file object (dataframe)
#' @return A nucleotide depth object
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
              tidyr::complete(pos = 1:max(pos), fill = list(n = 0))) %>%
    tidyr::unnest()
}


#' Normalise nucleotide depth using a normalization map
#'
#' @export
#' @param depths A nucleotide depth object
#' @param qc_norm A normalization map
#' @return A normalized depth object
normalise_depth <- function(depths, qc_norm)
{
  depths %>%
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
    dplyr::summarise(n = mean(n)) %>%
    dplyr::ungroup()
}


#' Compute moving average
#'
#' @export
#' @param depths A nucleotide depth object
#' @param window Half-size of the window
MA <- function(depths, window)
{
  depths %>%
    dplyr::filter(n > 0) %>%
    split(.$genotype) %>%
    lapply(function(x)
           {
             x %>%
               dplyr::mutate(roll = pos %>% purrr::map_dbl(~mean(x$n[x$pos %in% seq(. - window, . + window)]))) %>%
               tidyr::complete(pos = 1:max(pos), fill = list(genotype = x$genotype %>% unique,roll = 0, n = 0))
           }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(n = roll)
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
#' @export
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
#' @export
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


#' Plot extracted features
#'
#' @export
#' @param features A sam object (dataframe) with parsed CIGARs and extracted features
plot_features <- function(features)
{
  ggplot2::ggplot(features) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = feature_pos, color = feature, alpha = n)) +
    ggplot2::scale_alpha_continuous(range = 0:1) +
    ggplot2::scale_color_manual(values = c("left" = "red", "right" = "blue")) +
    ggplot2::scale_x_continuous(limits = c(0, 8000)) +
    ggplot2::facet_grid(~genotype) +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.text = ggplot2::element_blank())
}


#' Extract unaligned sequences
#'
#' @export
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


#' Clean a raw blat file
#'
#' @export
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
  dplyr::summarise(match = max(match)) %>%
  dplyr::ungroup()
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
    dplyr::select(-n)
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
                     dplyr::filter(feat == 2))
}


#' Tag a blat object with quality measures
#'
#' @export
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
    dplyr::filter(quality == max(quality)) %>%
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
#' @export
#' @param blat A tagged blat file
#' @return A summary table of results
summarise_blat <- function(blat)
{
  blat %>%
    dplyr::group_by(genotype, feature, position, chr, chr_position) %>%
    dplyr::add_count() %>%
    dplyr::summarise(n = max(n), quality = max(quality), match = max(match)) %>%
    dplyr::arrange(dplyr::desc(quality), dplyr::desc(n), dplyr::desc(match))
}
