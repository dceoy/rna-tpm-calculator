#!/usr/bin/env Rscript

'Tally and transform result data

Usage: create_matrix.R [--output <path>]

Options:
-h, --help      Print help and exit
--output        Pass an output directory path' -> doc

sapply(c('docopt', 'stringr', 'tidyverse'), require, character.only = TRUE)

main <- function(opts) {
  output_dir <- ifelse(is.null(opts[['--output']]),
                       './output/',
                       ifelse(str_detect(opts[['--out']], '/$'),
                              opts[['--output']],
                              str_c(opts[['--output']], '/')))

  df_count <- read_csv(file = str_c(output_dir, 'summary/read_count.csv')) %>%
    mutate(src_type = str_c(src, '_', type)) %>%
    select(id, src_type, count) %>%
    spread(src_type, count) %>%
    mutate(bam_unmapped_uniq = bam_total - bam_mapped,
           bam_uniq = bam_total - bam_secondary,
           bam_mapped_uniq = bam_mapped - bam_secondary) %>%
    select(id, bam_uniq, bam_mapped_uniq, bam_unmapped_uniq, bam_mapped, bam_secondary, bam_total, fastq_qc, fastq_raw)

  return(c(list(read_count_unformatted = write_csv(df_count,
                                                   path = str_c(output_dir, 'summary/read_count_matrix_unformatted.csv')),
                read_count_formatted = write_csv(select(df_count,
                                                        cell_id = id,
                                                        fastq = fastq_raw,
                                                        qc_passed = fastq_qc,
                                                        mapped = bam_mapped_uniq,
                                                        unmapped = bam_unmapped_uniq),
                                                 path = str_c(output_dir, 'summary/read_count_matrix.csv'))),
           lapply(list(gene_tpm = list(file = 'rsem_aligned.genes.results', tag = 'gene_id'),
                       transcript_tpm = list(file = 'rsem_aligned.isoforms.results', tag = 'transcript_id')),
                  function(l, sample_dir) {
                    write_csv(spread(bind_rows(lapply(list.files(sample_dir),
                                                      function(dn) {
                                                        return(mutate(select(read_tsv(file = str_c(sample_dir, dn, l$file, sep = '/')),
                                                                             matches(l['tag']), TPM),
                                                                      sample_id = dn))
                                                      })),
                                     sample_id, TPM, fill = NA),
                              path = str_c(output_dir, 'summary/tpm_', l['tag'], '.csv'))
                  },
                  sample_dir = str_c(output_dir, 'sample'))))
}

main(opts = docopt::docopt(doc))
