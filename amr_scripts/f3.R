

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

an_glass_df_pre <- rbind(an_glass_df_1,an_glass_df_2) %>%
  mutate(ab = trimws(as.character(ab))) %>%
  distinct()

glass_opts_join <- glass_opts %>%
  dplyr::select(antibiotic,antibiotic_class_glass) %>%
  distinct() %>%
  mutate(antibiotic = trimws(as.character(antibiotic)))

an_glass_df <- an_glass_df_pre %>%
  left_join(glass_opts_join, by = c('ab' = 'antibiotic')) %>%
  mutate(ab_class = antibiotic_class_glass) %>%
  distinct()



priority_glass_pathogens <-unique(an_glass_df$mo_organism)


message('Beginning GLASS analysis ....')
cat('Beginning GLASS analysis ....\n')




for (y in unique(an_df_long$yr)) {


for (i in par_df$id) {

  par=par_df$param[par_df$id==i]

  par_var_name=par_df$var_name[par_df$id==i]


  cntry = cntry

  par=par

  par_var_name=par_var_name

 #spaecimen distribution
  if (y==unique(an_df_long$yr)[1] & par==unique(par_df$param)[1]){

  specimen_distribution(an_glass_df)
  }else{
    NULL
  }


  abs_ref <- unique(an_glass_df$ab)


  #foreach (px=priority_glass_pathogens,.packages=loadedNamespaces(), .export = vars_fun) %dopar% {

  for(i in seq_along(priority_glass_pathogens)){
    org_name <- priority_glass_pathogens[i]
    org_name_dir <-str_replace_all(org_name," ","_")
    org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
    org_res_dir <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir, y,'Identified pathogen')
    org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,y, 'Identified pathogen',par)
    org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,'Trends','Identified pathogen', 'year')
    org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,'Trends','Identified pathogen',par, 'year')

    if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
    if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
    if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
    if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

    amr_individual_pathogens(an_glass_df,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                             org_res_dir_trends,org_res_dir_trends_par)
# }
  }

 # stopCluster(cl)


  org_name = "Staphylococcus aureus"
  org_name_dir <-paste(org_name,"mrsa")
  org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
  org_res_dir <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir, y,'Identified pathogen')
  org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,y, 'Identified pathogen',par)
  org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,'Trends','Identified pathogen', 'year')
  org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,'Trends','Identified pathogen',par, 'year')

  if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
  if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
  if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
  if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

  mrsa_analysis(an_glass_df,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                           org_res_dir_trends,org_res_dir_trends_par)


 # foreach (px=priority_glass_pathogens,.packages=loadedNamespaces(), .export = vars_fun) %dopar% {

  for(i in seq_along(priority_glass_pathogens)){
    org_name <- priority_glass_pathogens[i] ##extracting the genus name
    org_name_dir <-str_replace_all(org_name," ","_")
    org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
    org_res_dir <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,y,"pathogen_grp")
    org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,y, 'pathogen_grp',par)
    org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,'Trends','pathogen_grp', 'year')
    org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/GLASS"),org_name_dir,'Trends','pathogen_grp',par, 'year')

    if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
    if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
    if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
    if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

    amr_pathogen_groups(an_glass_df,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                        org_res_dir_trends, org_res_dir_trends_par)

  }
  # }
  #
  # stopCluster(cl)

}
  message(paste0(y),' done ....')
  cat(paste0(y),' done ....\n')
}


message('GLASS Analysis successfully completed ....')
cat(' GLASS Analysis successfully completed ....\n')
