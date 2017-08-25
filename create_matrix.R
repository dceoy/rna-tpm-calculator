#!/usr/bin/env Rscript

'Tally and transform result data

Usage: create_matrix.R [--output=<path>]

Options:
-h, --help        message help and exit
--output=<path>   Pass an output directory path' -> doc

main <- function(opts) {
  suppressMessages(sapply(c('stringr', 'tidyverse'), require, character.only = TRUE))
  extract_read_count_df <- function(file) {
    return(select(mutate(spread(select(mutate(read_csv(file = file),
                                              src_type = str_c(src, '_', type)),
                                       id, src_type, count),
                                src_type, count),
                         bam_unmapped_uniq = bam_total - bam_mapped,
                         bam_uniq = bam_total - bam_secondary,
                         bam_mapped_uniq = bam_mapped - bam_secondary),
                  id, bam_uniq, bam_mapped_uniq, bam_unmapped_uniq,
                  bam_mapped, bam_secondary, bam_total, fastq_qc, fastq_raw))
  }
  extract_tpm_df <- function(basefile, id_type, sample_dir) {
    spread(bind_rows(lapply(list.files(sample_dir),
                            function(dn) {
                              return(mutate(select(read_tsv(file = str_c(sample_dir, dn, basefile, sep = '/')),
                                                   matches(id_type), TPM),
                                            sample_id = dn))
                            })),
           sample_id, TPM, fill = NA)
  }

  output_dir <- ifelse(is.null(opts[['--output']]),
                       './output/',
                       ifelse(str_detect(opts[['--output']], '/$'),
                              opts[['--output']],
                              str_c(opts[['--output']], '/')))
  sample_dir <- str_c(output_dir, 'sample/')
  summary_dir <- str_c(output_dir, 'summary/')

  message('1. unformatted read count table')
  write_csv(suppressMessages(df_count <- extract_read_count_df(file = str_c(summary_dir, 'read_count.csv'))),
            path = str_c(summary_dir, 'read_count_matrix_unformatted.csv'))

  message('2. formatted read count table')
  write_csv(select(df_count,
                   sample_id = id,
                   fastq = fastq_raw,
                   qc_passed = fastq_qc,
                   mapped = bam_mapped_uniq,
                   unmapped = bam_unmapped_uniq),
            path = str_c(summary_dir, 'read_count_matrix.csv'))

  message('3. transcript tpm table')
  write_csv(suppressMessages(extract_tpm_df(id_type = 'transcript_id',
                                            basefile = 'rsem_aligned.isoforms.results',
                                            sample_dir = sample_dir)),
            path = str_c(summary_dir, 'tpm_transcript.csv'))

  message('4. gene tpm table')
  write_csv(suppressMessages(extract_tpm_df(id_type = 'gene_id',
                                            basefile = 'rsem_aligned.genes.results',
                                            sample_dir = sample_dir)),
            path = str_c(summary_dir, 'tpm_gene.csv'))

  message('done.')
}

suppressMessages(require('docopt'))
main(opts = docopt::docopt(doc))
