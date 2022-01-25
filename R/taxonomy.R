

#' Read lineages*.csv file and convert it to compressed R object
#'
#' This will be used for taxonomy matching.
lineagecsv_to_robj = function (lincsv, out) {
   taxonomy.dt = data.table::fread(lincsv, select=1:8)
   taxonomy.dt = taxonomy.dt[genus != '']
   taxonomy.dt[, species := gsub('^Candidatus ', '', species)]
   taxonomy.dt = unique(taxonomy.dt)
   saveRDS(taxonomy.dt, file=out, compress=T)
}

#' Load RExMap taxonomy reference file
#'
#' This is a NCBI Taxonomy database, exported as an compressed R object.
#'
#' @param tax_file Full path to the Rdata taxonomy file. Default:
#' rexmap_option('taxonomy_file').
#'
#' @export
load_taxonomy = function (tax_file=rexmap_option('taxonomy_file')) {
   readRDS(tax_file)
}

#' Parse genus names and strain counts from each OSU species string
#'
osuab_genuses = function (osuab, ws='_', split_char=',') {
  dt = osuab[, {
    species_list = strsplit(species[1], split_char, fixed=T)[[1]]
    genus_counts = lapply(species_list, function (s) {
      genus = gsub(paste0('^([^', ws, ']+)', ws, '.*'), '\\1', s)
      if (genus == 'Candidatus') {
         # For "Candidatus" genus, actual genus name is usually the second
         # word.
         genus = gsub(
            paste0('^[^', ws, ']+', ws, '([^', ws, ']+)', ws, '.*'),
            '\\1', s
         )
      }
      strain_count = gsub(paste0('.*', ws, '\\[([0-9]+)\\]$'), '\\1', s)
      if (strain_count == s) strain_count = '1'
      return(c(genus, strain_count))
    })
    genus_counts_unique = list()
    for (i in 1:length(genus_counts)) {
      genus_i = genus_counts[[i]][1]
      counts_i = as.integer(genus_counts[[i]][2])
      if (genus_i %in% names(genus_counts_unique)) {
        genus_counts_unique[[genus_i]] = genus_counts_unique[[genus_i]] + counts_i
      } else {
        genus_counts_unique[[genus_i]] = counts_i
      }
    }
    # Now pool together same genuses and add up the counts
    list('genus'=names(genus_counts_unique),
         # 'order_from_genus'=as.character(NA),
         'strain_count'=as.integer(unname(genus_counts_unique)),
         'pctsim'=unique(pctsim))
  }, by=.(osu_id)]
  # dt[genus %like% 'les$', c('order_from_genus', 'genus') := list(genus, NA)]
  return(dt)
}

#' Add NCBI taxonomy to an OSU table
#'
#' @param table Data table containing columns 'osu_id', 'pctsim' and 'species'.
#' Typically used with the output table from then \code{\link{abundance}}
#' function.
#'
#' @export
taxonomy = function (osu_table, verbose=rexmap_option('verbose'), show_count=TRUE,
                     ws='_', split_char=',') {
  if (show_count) {
    sep_str = '_['
    col_str = '],'
    fin_str = ']'
  } else {
    sep_str = ''
    col_str = ','
    fin_str = ''
  }

  if (verbose) cat('* load taxonomy...')
  taxonomy.dt = load_taxonomy()
  if (verbose) cat('OK.\n* generating unique ranks...')
  taxonomy.dt[, c('tax_id', 'species') := NULL]
  taxonomy.dt = unique(taxonomy.dt)
  tax_family.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class, order, family)])
  tax_order.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class, order, family=NA)])
  tax_class.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class, order=NA, family=NA)])
  tax_phylum.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class=NA, order=NA, family=NA)])
  if (verbose) cat('OK.\n* extracting genus/species from OSU table...')
  osu_abundance_table = unique(osu_table[, .(osu_id, pctsim, species)])
  osu_ab_g.dt = osuab_genuses(osu_abundance_table, ws=ws, split_char=split_char)
  if (verbose) cat('OK.\n* matching genus...')
  osu_ab_g2.dt = merge(osu_ab_g.dt, taxonomy.dt, by='genus', all.x=T) # Genus matches
  col_order = c('osu_id', 'strain_count', 'pctsim', 'superkingdom', 'phylum',
                'class', 'order', 'family', 'genus')
  setcolorder(osu_ab_g2.dt, col_order)
  if (verbose) cat('OK.\n* matching family...')
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, family],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_family.dt,
                            by.x='genus', by.y='family')
                 list(dt[, superkingdom], dt[, phylum], dt[, class], dt[, order],
                      dt[, genus], NA)
               }]
  if (verbose) cat('OK.\n* matching order...')
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, order],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_order.dt,
                            by.x='genus', by.y='order')
                 list(dt[, superkingdom], dt[, phylum], dt[, class], dt[, genus],
                      NA, NA)
               }]
  if (verbose) cat('OK.\n* matching class...')
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, class],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_class.dt,
                            by.x='genus', by.y='class')
                 list(dt[, superkingdom], dt[, phylum], dt[, genus],
                      NA, NA, NA)
               }]
  if (verbose) cat('OK.\n* matching phylum...')
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, phylum],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_phylum.dt,
                            by.x='genus', by.y='phylum')
                 list(dt[, superkingdom], dt[, genus],
                      NA, NA, NA, NA)
               }]
  # Now count unique taxonomic ranks for each osu_id
  if (verbose) cat('OK.\n* counting uniques...')
  if (show_count) {
     osu_ab_ranks.dt = osu_ab_g2.dt[, {
        kingdom_num = .SD[!is.na(superkingdom), .(kingdom_num = sum(strain_count)),
                          by=superkingdom][order(-kingdom_num)][
                             , paste(
                                superkingdom, kingdom_num,
                                sep=sep_str, collapse=col_str)
                             ]
        phylum_num = .SD[!is.na(phylum), .(phylum_num = sum(strain_count)),
                         by=phylum][order(-phylum_num)][
                            , paste(
                               phylum, phylum_num,
                               sep=sep_str, collapse=col_str)
                            ]
        class_num = .SD[!is.na(class), .(class_num = sum(strain_count)),
                        by=class][order(-class_num)][
                           , paste(
                              class, class_num,
                              sep=sep_str, collapse=col_str)
                           ]
        order_num = .SD[!is.na(order), .(order_num = sum(strain_count)),
                        by=order][order(-order_num)][
                           , paste(
                              order, order_num,
                              sep=sep_str, collapse=col_str)
                           ]
        family_num = .SD[!is.na(family), .(family_num = sum(strain_count)),
                         by=family][order(-family_num)][
                            , paste(
                               family, family_num,
                               sep=sep_str, collapse=col_str)
                            ]
        genus_num = .SD[!is.na(genus) & !grepl('^[a-z]', genus),
                        .(genus_num = sum(strain_count)),
                        by=genus][order(-genus_num)][
                           , paste(
                              genus, genus_num,
                              sep=sep_str, collapse=col_str)
                           ]
        list(
           'kingdom'=paste0(kingdom_num, fin_str),
           'phylum'=paste0(phylum_num, fin_str),
           'class'=paste0(class_num, fin_str),
           'order'=paste0(order_num, fin_str),
           'family'=paste0(family_num, fin_str),
           'genus'=paste0(genus_num, fin_str))
     }, by=.(osu_id, pctsim)]
  } else {
     osu_ab_g2.dt = osu_ab_g2.dt[order(osu_id, -pctsim, genus)]
     osu_ab_ranks.dt = osu_ab_g2.dt[
        ,
        list(
           'kingdom'=paste(unique(na.omit(superkingdom)), collapse=split_char),
           'phylum'=paste(unique(na.omit(phylum)), collapse=split_char),
           'class'=paste(unique(na.omit(class)), collapse=split_char),
           'order'=paste(unique(na.omit(order)), collapse=split_char),
           'family'=paste(unique(na.omit(family)), collapse=split_char),
           'genus'=paste(unique(na.omit(genus)), collapse=split_char)
        )
        , by=.(osu_id, pctsim)
     ]
  }
  osu_ab_ranks.dt[!(grepl('^[A-Z]', kingdom)), kingdom := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', phylum)), phylum := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', class)), class := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', order)), order := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', family)), family := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', genus)), genus := NA]
  if (verbose) cat('OK.\n* Done.\n')
  return(osu_ab_ranks.dt[])
}

