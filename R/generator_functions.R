#' Generate a global_align function
#'
#' @export
#' @param HPVs HPVs reference file
#' @param fastqdir Directory of fastq files
#' @param threads Number of threads to use for bowtie2
#' @param write Write the sam file on disk?
#' @param dir Directory in which to write the global sam file if write == T
#' @return A function accepting a fastq filename prefix that returns a parsed sam object (dataframe) for global alignment
Gglobal_align <- function(HPVs, fastqdir, threads = NA, write = F, dir = NULL)
{
  if (write)
    dir.create(paste0(fastqdir, "/", dir), showWarnings = F, recursive = T)

  function(fastq)
  {
    bowtie_global_align(HPVs, paste0(fastqdir, "/", fastq), threads) -> sam

    if (write)
      cat(sam, file = paste0(fastqdir, "/", dir, "/", fastq, ".sam"), sep = "\n")

    sam %>%
      parse_sam %>%
      dplyr::select(read = qname, genotype = rname, pos, cigar, seq)
  }
}


#' Generate a local_align function
#'
#' @export
#' @param HPVdir A directory of HPV reference files
#' @param fastqdir Directory of fastq files
#' @param threads Number of threads to use for bowtie2
#' @param write Write the sam file on disk?
#' @param dir Directory in which to write the local sam file if write == T
#' @return A function accepting a genotype name and fastq filename prefix that returns a parsed sam object (dataframe) for local alignment
Glocal_align <- function(HPVdir, fastqdir, threads = NA, write = F, dir = NULL)
{
  if (write)
    dir.create(paste0(fastqdir, "/", dir), showWarnings = F, recursive = T)

  function(genotype, fastq)
  {
    bowtie_local_align(paste0(HPVdir, "/", genotype), paste0(fastqdir, "/", fastq), threads) -> sam

    if (write)
      cat(sam, file = paste0(fastqdir, "/", dir, "/", fastq, "-", genotype, ".sam"), sep = "\n")

    sam %>%
      parse_sam %>%
      dplyr::select(read = qname, pos, cigar, seq)
  }
}


#' Select genotype
#'
#' @export
#' @param ... Filtering conditions to keep a genotype
#' @return A function accepting a sam object for global alignment and returns a list of selected genotypes
Gselect_genotype <- function(...)
{
  function(globalsam)
  {
    globalsam %>%
      dplyr::count(genotype, sort = T) %>%
      tidyr::separate(genotype, c("genotype", "variant"), sep = "_") %>%
      dplyr::group_by(genotype) %>%
      dplyr::top_n(1, n) %>%
      dplyr::distinct(genotype, .keep_all = T) %>%
      tidyr::unite(genotype, genotype, variant) %>%
      dplyr::filter(...) %>%
      dplyr::mutate(genotype = genotype %>% sub(pattern = "_NA", replacement = "")) %>%
      dplyr::select(-n)
  }
}


#' Generate a plot function
#'
#' @export
#' @param write Write plot on disk?
#' @param fastqdir Directory of fastq files
#' @param dir Directory in which to write the plot if write == T
#' @return A function accepting a sam object and fastq name and returns a plot object
Gplot_reads <- function(write = F, fastqdir = NULL, dir = NULL)
{
  if (write)
    dir.create(paste0(fastqdir, "/", dir), showWarnings = F, recursive = T)

  function(sam, fastq = NULL)
  {
    ggplot_reads(sam) -> Plot

    if (write)
      ggplot2::ggsave(plot = Plot, filename = paste0(fastqdir, "/", dir, "/", fastq, ".png"))

    Plot
  }
}


#' Generate a gap plot function
#'
#' @export
#' @param write Write plot on disk?
#' @param fastqdir Directory of fastq files
#' @param dir Directory in which to write the plot if write == T
#' @return A function accepting a local plot object, a features plot object,and a fastq name, and returning a plot object
Gplot_gaps <- function(write = F, fastqdir = NULL, dir = NULL)
{
  if (write)
    dir.create(paste0(fastqdir, "/", dir), showWarnings = F, recursive = T)

  function(local_plot, features_plot, fastq = NULL)
  {
    patchwork::wrap_plots(local_plot, features_plot, ncol = 1, heights = c(.9, .1)) -> Plot

    if (write)
      ggplot2::ggsave(plot = Plot, filename = paste0(fastqdir, "/", dir, "/", fastq, ".png"))
  }
}


#' Generate the write fasta function
#'
#' @export
#' @param fastqdir Directory of fastq files
#' @param dir Directory in which to write the fasta file
#' @return A function accepting a sam file with parsed CIGARs and extracted sequences, and a fastq name, which writes the fasta file on disk
Gwrite_fasta <- function(fastqdir, dir)
{
  dir.create(paste0(fastqdir, "/", dir), showWarnings = F, recursive = T)

  function(reads, fastq)
  {
    paste(paste0(">", reads$read, "|", reads$genotype, "|", reads$feature, "|", reads$feature_pos), reads$nalign_seq, sep = "\n", collapse = "\n") %>%
      cat(file = paste0(fastqdir, "/", dir, "/", fastq, ".fa"))
  }
}


#' Generate the blat function
#'
#' @export
#' @param fastqdir Directory of fastq files
#' @param arguments Additionnal arguments for command line blat
#' @param fasta_dir Directory of fasta files
#' @param blat_dir Directory in which to write the results from blat
#' @return A function accepting a fastq name, which runs blat on the corresponding fasta file and produces the resulting output from blat
Gblat <- function(fastqdir, arguments = NULL, fasta_dir, blat_dir)
{
  dir.create(paste0(fastqdir, "/", blat_dir), showWarnings = F, recursive = T)

  function(fastq)
  {
    fasta <- paste0(fastqdir, "/", fasta_dir, "/", fastq, ".fa")
    psl <- paste0(fastqdir, "/", blat_dir, "/", fastq, ".tsv")

    system2("blat",
            c("hg19.2bit",
              fasta,
              psl,
              "-noHead",
              arguments))
  }
}


#' Generate the read_blat function
#'
#' @export
#' @param fastqdir Directory of fastq files
#' @param blat_dir Directory in which to read the blat file
#' @return A function accepting a fastq name, which reads the corresponding blat file
Gread_blat <- function(fastqdir, blat_dir)
{
  function(fastq)
  {
    paste0(fastqdir, "/", blat_dir, "/", fastq, ".tsv") %>%
      utils::read.delim(col.names = c("match","mis-match","rep.match","N's","Q gap count","Q gap bases","T gap count","T gap bases","strand","Q name","Q size","Q start","Q end","T name","T size","T start","T end","block count","blockSizes","qStarts","tStarts"),
                 stringsAsFactors = F,
                 check.names = F,
                 header = F) %>%
      tidyr::separate(`Q name`, c("read", "genotype", "feature", "position"), sep = "\\|")
  }
}


#' Generate the write_blat function
#'
#' @export
#' @param fastqdir Directory of fastq files
#' @param blat_dir Directory in which to write the blat file
#' @return A function accepting a blat object (raw, cleaned, or summarised), a fastq name, and a suffix name, which writes the corresponding blat file with the suffix
Gwrite_blat <- function(fastqdir, blat_dir)
{
  function(blat, fastq, suffix)
  {
    blat %>%
      utils::write.table(paste0(fastqdir, "/", blat_dir, "/", fastq, "_", suffix, ".tsv"),
                  row.names = F,
                  sep = "\t",
                  quote = F)
  }
}
