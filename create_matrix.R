#!/usr/bin/env Rscript

'Tally and transform result data

Usage: create_matrix.R <type>

Options:
  -h, --help      Print help and exit
  <type>          Pass a matrix type {read,tpm}' -> doc

sapply(c('docopt', 'stringr', 'tidyverse'), require, character.only = TRUE)

main <- function(opts) {
  count_reads <- function() {
    df_count <- read_csv('summary/read_count.csv') %>%
      mutate(src_type = str_c(src, '_', type)) %>%
      select(id, src_type, count) %>%
      spread(src_type, count) %>%
      mutate(bam_unmapped_uniq = bam_total - bam_mapped,
             bam_uniq = bam_total - bam_secondary,
             bam_mapped_uniq = bam_mapped - bam_secondary) %>%
      select(id, bam_uniq, bam_mapped_uniq, bam_unmapped_uniq, bam_mapped, bam_secondary, bam_total, fastq_qc, fastq_raw)
    df_count %>%
      write_csv(file = 'summary/read_count_matrix_unformatted.csv', row.names = FALSE)
    df_count %>%
      select(cell_id = id,
             fastq = fastq_raw,
             qc_passed = fastq_qc,
             mapped = bam_mapped_uniq,
             unmapped = bam_unmapped_uniq) %>%
      write_csv(file = 'summary/read_count_matrix.csv', row.names = FALSE)
  }
  tally_tpm <- function() {
    lapply(list(c(path = 'rsem/GRCm38_p4_TCF7.genes.results', tag = 'gene_id'),
                c(path = 'rsem/GRCm38_p4_TCF7.isoforms.results', tag = 'transcript_id')),
           function(l, map_dir, tpm_dir) {
             write.csv(spread(bind_rows(lapply(str_replace(list.dirs(map_dir, recursive = FALSE),
                                                           str_c(map_dir, '/'),
                                                           ''),
                                               function(dn) {
                                                 str_c(map_dir, dn, l['path'], sep = '/') %>%
                                                   fread() %>%
                                                   select(matches(l['tag']), TPM) %>%
                                                   mutate(cell_id = dn) %>%
                                                   return()
                                               })),
                              cell_id, TPM, fill = NA),
                       row.names = FALSE, file = str_c(tpm_dir, '/tpm_', l['tag'], '.csv'))
           },
           map_dir = 'map',
           tpm_dir = tpm_dir)
  }
  return(switch(opts[['<type>']],
                'read' = count_reads(),
                'tpm' = tally_tpm()))
}

main(opts = docopt::docopt(doc))
