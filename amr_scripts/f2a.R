##Update Patient location information
file_path <- file.path(amr_updates_dir,"patient_location_type.xlsx")

if (file.exists(file_path)) {
 patient_loc_updates <- read.xlsx(file_path)%>% filter(!is.na(options))  # or read.csv, fread, etc.
} else {
  patient_loc_updates=loc_options%>% filter(!is.na(options))
  message("File not found — matching_GLASS_specimen_types.xlsx")
}

#
name_loc <- lkp_facility %>% left_join(patient_loc_updates %>% rename(location_type=options),
                           by=c('Patient Location Type'='my_dataset')) %>%
  dplyr::select(r_id,location_type)

# Start the downstream analysis -------------------------------------------

message('Beginning downstream analysis ....')
cat('Beginning downstream analysis ....\n')

# Please note Rates are only shown when n >=30 according to GLASS reccommendations

#create an analysis dataframe
#for glass options
spec_options <- read_excel('amr_resources/list_glass_2.xlsx', sheet = 'whonet_spec_codes')



an_df <- sir_outcomes_df_wide %>%
  left_join(lkp_demographics %>% dplyr::select(Age, Sex,r_id), by='r_id') %>%
  left_join(lkp_facility %>% dplyr::select(`Laboratory Name`,r_id), by='r_id') %>%
  ##cleanung Age cols
  mutate(Age_s=toupper(gsub("[[:digit:]]", "", Age)),
         Age_n= as.numeric(gsub("[^0-9.-]", "", Age)),
         Age_conv=ifelse(Age_s=='D',round(Age_n/365,0),
                         ifelse(Age_s=='W',round(Age_n/7,0),
                                ifelse(Age_s=='M',round(Age_n/12,0),
                                       ifelse(Age_s=='Y',round(Age_n/1,0),Age)))),
         Age_conv=as.numeric(Age_conv),
         Age_g=as.character(age_groups(Age_conv, split_at = c(1,5,20,50,65))),
         Age_g=ifelse(Age_g=='0', '<1',Age_g),
         Age_g=factor(Age_g, levels=c('<1','1-4','5-19', '20-49','50-64','65+')),
         Sex =str_to_title(Sex),
         specimen_type=tolower(specimen_type)) %>%
  left_join(spec_options, by=c('specimen_type'='Spec_Type_code')) %>%
  mutate(specimen_type=ifelse(!is.na(Spec_name), Spec_name, specimen_type)) %>%
  distinct(uid, `Laboratory Name`, specimen_date, specimen_type, organism, .keep_all = T)


##Saving the interpreted file
openxlsx::write.xlsx(an_df,file = paste0(cntry,"/Results_AMR/Interpreted.results.",date_var,".xlsx"))


ab_cols <- as.character(unique(sir_outcomes_df$ab))  #antibiotic columns

if(file.exists('amr_resources/ab_class_list.csv')){
  ab_class_list <- read.csv('amr_resources/ab_class_list.csv')
}else{
  stop("Error: The file ab_class_list.csv does not exist. Please add file to test-data and try again.")
}

an_df_long <- an_df %>%
  pivot_longer(cols=ab_cols, names_to = 'ab', values_to = 'interpreted_res') %>%
  filter(!is.na(interpreted_res)) %>%
  mutate(R=ifelse(interpreted_res=='R',1,0),       ##resistance column
         genus=str_split_i(mo_organism, ' ',1),
         yr_month=format(specimen_date, "%Y-%m"),
         yr=format(specimen_date, "%Y"),
         quarter= ifelse(grepl('-01|-02|-03', yr_month), 'Q1',
                         ifelse(grepl('-04|-05|-06', yr_month), 'Q2',
                                ifelse(grepl('-07|-08|-09', yr_month), 'Q3',
                                       ifelse(grepl('-10|-11|-12', yr_month), 'Q4',
                                              '')))),
         yr_quarter=paste(yr,quarter)) %>%
  mutate(interpreted_res=as.sir(interpreted_res, clean = TRUE)) %>%
  left_join(ab_class_list, by='ab', relationship = "many-to-many") %>%
  dplyr::select(-ab_selector, -X) %>%
  mutate(ab_class=ifelse(is.na(ab_class), ab_group(ab), ab_class)) %>%   #if class is not in the current list, source it
  distinct() %>% left_join(name_loc, by='r_id')


# Write data to file ------------------------------------------------------

openxlsx::write.xlsx(lkp_organisms,file = paste0(cntry,"/Results_AMR/Organisms.",date_var,".xlsx"))
openxlsx::write.xlsx(lkp_demographics,file = paste0(cntry,"/Results_AMR/Demographics.",date_var,".xlsx"))
openxlsx::write.xlsx(lkp_facility,file = paste0(cntry,"/Results_AMR/Facilities.",date_var,".xlsx"))
openxlsx::write.xlsx(sir_outcomes_df_wide,file = paste0(cntry,"/Results_AMR/AST.results.",date_var,".xlsx"))
openxlsx::write.xlsx(excluded_rec,file = paste0(cntry,"/Results_AMR/Intrinsic.noguidelines.results.",date_var,".xlsx"))


#openxlsx::write.xlsx(lkp_organisms,file = paste0("Organisms.",date_var,".xlsx"))
#openxlsx::write.xlsx(lkp_demographics,file = paste0("Demographics.",date_var,".xlsx"))
#openxlsx::write.xlsx(lkp_facility,file = paste0("Facilities.",date_var,".xlsx"))
#openxlsx::write.xlsx(sir_outcomes_df_wide,file = paste0("AST.results.",date_var,".xlsx"))
#openxlsx::write.xlsx(excluded_rec,file = paste0("Intrinsic.noguidelines.results.",date_var,".xlsx"))



# Get bug-drug combinations  ----------------------------------------------
# where at least 30 (default) isolates are available per species

#-Gilbeert: commenting this one out
#bug_drug_combos <- format(AMR::bug_drug_combinations(an_df_long))
#openxlsx::write.xlsx(bug_drug_combos,file = file.path("Results",paste0("Bug_drug_combinations.results.",date_var,".xlsx")))

# National antibiogram ----------------------------------------------------


#preprocessing to work in macs also

prep_for_antibiogram_mac <- function(df) {
  df %>%
    mutate(across(where(is_sir_eligible),
                  ~ suppressWarnings(as.sir(as.character(.)))))
}


# Detailed analysis -------------------------------------------------------

message('Profile analysis ....')
cat('Profile analysis ....\n')

## Define parameters for subgroup analyses
par_df <- tibble(param=c('Age', 'Sex', 'Specimen Type', 'location_type'#,
                        # 'Laboratory Name'
                         ), var_name=c('Age_g', 'Sex', 'specimen_type',
                                      'location_type'#, 'Laboratory Name'
                                      )) %>%
  mutate(id=paste0(param,var_name))

orgs_vec <- lkp_organisms %>% dplyr::pull(fullname)


if (exists("priority_bact_pathogens")) {
  priority_bact_vec <- orgs_vec[orgs_vec %in% priority_bact_pathogens]
} else {
  print("priority_bact_pathogens vector does NOT exist. Using universal priority pathogens list")

  priority_bact_pathogens <- c(
    "Acinetobacter baumannii",
    "Enterococcus faecalis",
    "Enterococcus faecium",
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Neisseria gonorrhoeae",
    "Pseudomonas aeruginosa",
    "Staphylococcus aureus",
    "Streptococcus agalactiae",
    "Streptococcus pneumoniae",
    "Streptococcus pyogenes",
    "Helicobacter pylori",
    "Mycobacterium tuberculosis"
  )

  priority_bact_pathogens_grp <- c("Citrobacter",
                                   "Enterobacter","Morganella","Salmonella",
                                   "Proteus","Serratia",
                                   "Shigella","Campylobacter")


  priority_bact_vec <- orgs_vec[orgs_vec %in% priority_bact_pathogens]

}

priority_fungi_pathogens   <-c( "Cryptococcus neoformans",  "Candida auris",
    "Aspergillus fumigatus",    "Candida albicans",   "Nakaseomyces glabrata",
    "Candida tropicalis",    "Candida parapsilosis",  "Pneumocystis jirovecii",
    "Talaromyces marneffei",  "Lomentospora prolificans")


priority_fungi_pathogens_grp <- c("Histoplasma","Paracoccidioides",
                                  "Scedosporium","Cladosporium",
                                  "Eumycetoma", "Mucorales",
                                  "Fusarium","Coccidioides")

#priority_fungal_vec <- priority_Fungi_pathogens


# Set up cluster

#registerDoParallel(cl, cores = detectCores() - 1)

#cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1), type = "PSOCK")
#doParallel::registerDoParallel(cl)

vars_fun <- c("an_df_long","y","par","par_df","an_df",
              #"an_glass_df",
              "indiv_ab_resistance", "overall_ab_resistance",
              "indiv_ab_resistance_sau","overall_ab_resistance_sau",
              "indiv_ab_resistance_genus", "overall_ab_resistance_genus",
              "antibiotic_classes_res_indiv", "antibiotic_classes_res_grp",
              "overall_resistance_plot","indiv_ab_resistance_plot",
              #"amr_individual_pathogens",
              "mrsa_analysis",
              #"amr_pathogen_groups",
              "indiv_ab_resistance_plot_trends", "overall_resistance_plot_trends",
              "indiv_ab_resistance_plot_trends_quarter", "overall_resistance_plot_trends_quarter"
)

# on.exit({
#   doParallel::registerDoSEQ()
#   try(parallel::stopCluster(cl), silent = TRUE)
# }, add = TRUE)
#
# # Dynamically block any connection/cluster objects from auto-export
# all_names <- ls(envir = .GlobalEnv, all.names = TRUE)
# connish <- Filter(function(nm) {
#   obj <- try(get(nm, envir = .GlobalEnv), silent = TRUE)
#   if (inherits(obj, "try-error")) return(FALSE)
#   inherits(obj, "connection") || inherits(obj, "cluster") || inherits(obj, "sockconn")
# }, all_names)


#doParallel::registerDoParallel(cl)
#foreach (y=an_df_long$yr,.packages=loadedNamespaces(), .export = vars_fun,.noexport = unique(c(connish, "cl", "con", "conn", "db", "drv")) , .errorhandling = "pass") %dopar% {

for (y in unique(an_df_long$yr)) {

for (i in par_df$id) {

  par=par_df$param[par_df$id==i]

  par_var_name=par_df$var_name[par_df$id==i]


  cntry = cntry

  par=par

  par_var_name=par_var_name

   # org_name='Escherichia coli',
  # abs_ref <- c("AMP", "CFR", "CTX", "CAZ", "CIP", "GEN", "TOB", "MEC", "MEM", "NIT", "TZP", "TMP", "SXT")


  abs_ref <- unique(an_df_long$ab)


  #Antibiogram

    an_dfx <- an_df %>% mutate(year=format(specimen_date, "%Y")) %>%
    filter(year==y)


  an_dfx <- prep_for_antibiogram_mac(an_dfx)


  ## ---- call per OS ----
  if (user_os == "Mac") {
    abg_df <- an_dfx %>%
      mutate_if(is_sir_eligible, as.sir) %>%
      mutate(bacteria = as.mo(bacteria)) %>%
      safe_antibiogram(antimicrobials = ab_cols)
  } else {
    abg_df <- an_dfx %>%
      mutate_if(is_sir_eligible, as.sir) %>%
      mutate(bacteria = as.mo(bacteria)) %>%
      safe_antibiogram()
  }

  # Creating a workbook and worksheet
  wb <- createWorkbook()
  addWorksheet(wb, "antibiogram")

  # Write data to Excel
  writeData(wb, sheet = "antibiogram", x = abg_df)

  # Apply color scale to 'Score' column (column B)
  conditionalFormatting(
    wb,
    sheet = "antibiogram",
    cols = 2:length(names(abg_df)),
    rows = 2:(nrow(abg_df)+1),  # +1 to account for header
    style = c("red", "yellow", "green"),
    type = "colourScale"
  )

  # Save to file
  saveWorkbook(wb, file.path(paste0(cntry,"/Results_AMR/National.antibiogram.results.",y,'_',date_var,".xlsx")), overwrite = TRUE)




  #Bacterial segment

  #doParallel::registerDoParallel(cl)
  #registerDoParallel(cl, cores = detectCores() - 1)
  #foreach (px=priority_bact_vec,.packages=loadedNamespaces(), .export = vars_fun, .errorhandling = "pass") %dopar% {

     for(i in seq_along(priority_bact_pathogens)){
    org_name <- priority_bact_pathogens[i]
    org_name_dir <-str_replace_all(org_name," ","_")
    org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
    org_res_dir <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir, y,'Identified pathogen')
    org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,y, 'Identified pathogen',par)
    org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,'Trends','Identified pathogen', 'year')
    org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,'Trends','Identified pathogen',par, 'year')

    if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
    if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
    if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
    if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

    amr_individual_pathogens(an_df_long,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                             org_res_dir_trends,org_res_dir_trends_par)
     }
 # }
  #stopCluster(cl)



  org_name = "Staphylococcus aureus"
  org_name_dir <-paste(org_name,"mrsa")
  org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
  org_res_dir <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir, y,'Identified pathogen')
  org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,y, 'Identified pathogen',par)
  org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,'Trends','Identified pathogen', 'year')
  org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,'Trends','Identified pathogen',par, 'year')

  if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
  if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
  if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
  if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

  mrsa_analysis(an_df_long,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                org_res_dir_trends,org_res_dir_trends_par)




 # registerDoParallel(cl, cores = detectCores() - 1)
  #foreach (px=priority_bact_vec,.packages=loadedNamespaces(), .export = vars_fun) %dopar% {

      for(i in seq_along(priority_bact_pathogens_grp)){
    org_name <- priority_bact_pathogens_grp[i] ##extracting the genus name
    org_name_dir <-str_replace_all(org_name," ","_")
    org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
    org_res_dir <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,y,"pathogen_grp")
    org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,y, 'pathogen_grp',par)
    org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,'Trends','pathogen_grp', 'year')
    org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/Bacteria"),org_name_dir,'Trends','pathogen_grp',par, 'year')

    if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
    if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
    if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
    if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

    amr_pathogen_groups(an_df_long,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                        org_res_dir_trends, org_res_dir_trends_par)

    #  }
  }

 # stopCluster(cl)



  #Fungal#priority_fungi_vec

  # registerDoParallel(cl, cores = detectCores() - 1)
  # foreach (px=priority_fungal_vec,.packages=loadedNamespaces(), .export = vars_fun) %dopar% {
  #
     for(i in seq_along(priority_fungi_pathogens)){
    org_name <- priority_fungi_pathogens[i]
    org_name_dir <-str_replace_all(org_name," ","_")
    org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
    org_res_dir <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir, y,'Identified pathogen')
    org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,y, 'Identified pathogen',par)
    org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,'Trends','Identified pathogen', 'year')
    org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,'Trends','Identified pathogen',par, 'year')

    if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
    if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
    if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
    if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

    amr_individual_pathogens(an_df_long,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                             org_res_dir_trends,org_res_dir_trends_par)
     }
  # }
  # stopCluster(cl)




 # registerDoParallel(cl, cores = detectCores() - 1)
  #foreach (px=priority_fungal_vec,.packages=loadedNamespaces(), .export = vars_fun) %dopar% {

      for(i in seq_along(priority_fungi_pathogens_grp)){
    org_name <- priority_fungi_pathogens_grp[i] ##extracting the genus name
    org_name_dir <-str_replace_all(org_name," ","_")
    org_name_dir <-str_replace_all(org_name_dir,"\\(|\\)","_")
    org_res_dir <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,y,"pathogen_grp")
    org_res_dir_par <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,y, 'pathogen_grp',par)
    org_res_dir_trends <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,'Trends','pathogen_grp', 'year')
    org_res_dir_trends_par <- file.path(paste0(cntry,"/Results_AMR/Fungi"),org_name_dir,'Trends','pathogen_grp',par, 'year')

    if(!dir.exists(org_res_dir)){dir.create(org_res_dir, recursive = T)}
    if(!dir.exists(org_res_dir_par)){dir.create(org_res_dir_par, recursive = T)}
    if(!dir.exists(org_res_dir_trends)){dir.create(org_res_dir_trends, recursive = T)}
    if(!dir.exists(org_res_dir_trends_par)){dir.create(org_res_dir_trends_par, recursive = T)}

    amr_pathogen_groups(an_df_long,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name,
                        org_res_dir_trends, org_res_dir_trends_par)

    #  }
  }


  }

  message(paste0(y),' done ....')
  cat(paste0(y),' done ....\n')
}
#stopCluster(cl)

##GLASS calculations
glass_opts <- read_excel('amr_resources/list_glass_2.xlsx') %>%
  mutate(mo_organism=as.mo(pathogen),
         antibiotic=as.ab(antibiotic),
         rec_combo=paste0(mo_organism,specimen,antibiotic))

##remap the specimen types in the analysis df

#unique specimen names
specimen_check_df <- data.frame(my_dataset=unique(an_df$specimen_type),
                                matched=c(rep('', length(unique(an_df$specimen_type))))) %>%
  filter(!(my_dataset %in% glass_opts$specimen))



message('Analysis successfully completed ....')
cat('Analysis successfully completed ....\n')
