

#After matching
#read in the user-validated file
file_path <- file.path(amr_updates_dir,"matching_GLASS_specimen_types.xlsx")

if (file.exists(file_path)) {
  glass_spec_updates <- read.xlsx(file_path)%>% filter(!is.na(matched))  # or read.csv, fread, etc.
} else {
  glass_spec_updates=specimen_check_df%>% filter(!is.na(matched))
  message("File not found — matching_GLASS_specimen_types.xlsx")
}

#updating the specimen names
for (i in glass_spec_updates$my_dataset) {
  an_df_long$specimen_type[an_df_long$specimen_type==i]=glass_spec_updates$matched[glass_spec_updates$my_dataset==i]

}

an_glass_df_1 <- an_df_long %>%  mutate(rec_combo=paste0(bacteria,specimen_type,ab)) %>%
  filter(rec_combo %in% glass_opts$rec_combo)

an_glass_df_2 <- an_df_long %>%  mutate(rec_combo=paste0(as.mo(genus),specimen_type,ab)) %>%  #for acinetobacter, salmonella, shigella
  filter(rec_combo %in% glass_opts$rec_combo) %>%
  mutate(mo_organism=genus)

an_glass_df <- rbind(an_glass_df_1,an_glass_df_2) %>%
  left_join(glass_opts %>% dplyr::select(antibiotic,antibiotic_class_glass ) %>% distinct() %>% mutate(antibiotic=as.character(antibiotic)),
            by=c('ab'='antibiotic')) %>%
  mutate(ab_class=antibiotic_class_glass) %>%
  distinct()



priority_glass_pathogens <-unique(an_glass_df$mo_organism)


message('Beginning GLASS analysis ....')
cat('Beginning GLASS analysis ....\n')


# route worker stdout/stderr to the master (prevents dangling/sunk connections)
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1),
                            outfile = NULL,                           # don't write to stdout/files
                            type = "PSOCK",
                            rscript_args = c("--vanilla"))             # clean worker sessions)


  ## ---- loops ----
  for (y in unique(data_yrs$yr)) {

    for (i in par_df$id) {

      par <- par_df$param[par_df$id == i]
      par_var_name <- par_df$var_name[par_df$id == i]

      #spaecimen distribution
      if (y==unique(an_df_long$yr)[1] & par==unique(par_df$param)[1]){

        specimen_distribution(an_glass_df)
      }else{
        NULL
      }

      abs_ref <- unique(an_glass_df$ab)


      message(paste0(y), " ",par," analysis now in progress...please hang in there ....")
      cat(paste0(y), " ",par," analysis now in progress...please hang in there ....\n")

      ## ---- Bacterial: individual pathogens ----
      # (Keeping re-registration; it’s okay and idempotent.)
      doParallel::registerDoParallel(cl)
      parallel::clusterCall(cl, function() options(warn = -1))

      foreach(
        px = priority_glass_pathogens,
        .packages = pkg_slim,
        .export   = unique(c(vars_fun, "y","par","par_var_name","abs_ref","cntry")),
        .noexport = unique(c(connish, "cl","con","conn","db","drv")),
        .errorhandling = "pass"
      ) %dopar% {

        parallel_guard()

        # fully qualify stringr calls so workers don’t depend on search path
        org_name <- px
        org_name_dir <- stringr::str_replace_all(org_name, " ", "_")
        org_name_dir <- stringr::str_replace_all(org_name_dir, "\\(|\\)", "_")

        org_res_dir         <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, y, "Identified pathogen")
        org_res_dir_par     <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, y, "Identified pathogen", par)
        org_res_dir_trends  <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, "Trends", "Identified pathogen", "year")
        org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, "Trends", "Identified pathogen", par, "year")

        mkpath(org_res_dir)
        mkpath(org_res_dir_par)
        mkpath(org_res_dir_trends)
        mkpath(org_res_dir_trends_par)

        amr_individual_pathogens(
          an_glass_df, org_res_dir, org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
          org_res_dir_trends, org_res_dir_trends_par
        )
        NULL  #  return nothing to avoid serializing large objects
      }

      ## ---- MRSA ----
      org_name <- "Staphylococcus aureus"
      org_name_dir <- paste(org_name, "mrsa")
      org_name_dir <- stringr::str_replace_all(org_name_dir, "\\(|\\)", "_")
      org_res_dir <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, y, "Identified pathogen")
      org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, y, "Identified pathogen", par)
      org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, "Trends", "Identified pathogen", "year")
      org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"), org_name_dir, "Trends", "Identified pathogen", par, "year")

      mkpath(org_res_dir)
      mkpath(org_res_dir_par)
      mkpath(org_res_dir_trends)
      mkpath(org_res_dir_trends_par)

      mrsa_analysis(
        an_glass_df, org_res_dir, org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
        org_res_dir_trends, org_res_dir_trends_par
      )

    } # end par loop

    message(paste0(y), " done ....")
    cat(paste0(y), " done ....\n")
  } # end year loop



message('GLASS Analysis successfully completed ....')
cat(' GLASS Analysis successfully completed ....\n')
