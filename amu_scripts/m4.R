

cat("Initializing updates on the dataset\n")
message("Initializing updates on the dataset\n")

##DDD updates
#read in the user-validated file
file_path <- file.path(amu_updates_dir,'DDD_information_updates.xlsx')

if (file.exists(file_path)) {
  ddd_updates <- read.xlsx(file_path)%>% rename(`atc_level_name`=atc_level_name)  # or read.csv, fread, etc.
} else {
  ddd_updates=ddd_updates
  message("File not found — ddd_info_updates")
}

#options to look up
ddd_updates=ddd_updates %>% mutate(name_route=paste0(`atc_level_name`, '_',tolower(Adm.R)))


ddd_ref_update <- bind_rows_match_classes(list(ddd_ref,ddd_updates %>% filter(DDD!=' ')))

#updating the atcs cleaned file
writexl::write_xlsx(ddd_ref_update, 'amc_resources/ab_molecules_amc.xlsx')


##Strength units
#read in the user-validated file
file_path <- file.path(amu_updates_dir,'strength_units_updates_b.xlsx')

if (file.exists(file_path)) {
  units_updates <- read.xlsx(file_path)  # or read.csv, fread, etc.
} else {
  units_updates=units_updates
  message("File not found — ddd_info_updates")
}

units_ref_update <- rbind(units_ref,
                          units_updates %>% filter(standardized_units!='') %>%
                            mutate_all(as.character())) %>%
  mutate(strength_unit_in_dataset=tolower(strength_unit_in_dataset)) %>%
  distinct(strength_unit_in_dataset, .keep_all = T)

#updating the atcs cleaned file
writexl::write_xlsx(units_ref_update, 'amc_resources/strength_units_updates_b.xlsx')


#incorporating the added ddd info
amu_dataset_comb2 <- amu_dataset_comb_r %>%  filter(name_route %in% ddd_ref$name_route)

unused_records_ddd <- subset(amu_r1, !(amu_r1$uid %in%  c(unique(amu_dataset1$uid),unique(amu_dataset_inhibitors$uid),
                                                      unique(amu_dataset_comb1$uid), unique(amu_dataset_comb2$uid))))

unused_records_amu <- subset(amu_r1, !(amu_r1$uid %in%  c(unique(amu_dataset1$uid),unique(amu_dataset_inhibitors$uid),
                                                          unique(amu_dataset_comb1$uid), unique(amu_dataset_comb2$uid),
                                                          unique(amu_dataset_comb_r_valid$uid))))

unused_records_without_date <- amu_raw_meta %>% filter(is.na(date_of_data_collection))



writexl::write_xlsx(unused_records_ddd, paste0(amu_updates_dir,'/unused_data_for_ddd.xlsx'))
writexl::write_xlsx(unused_records_amu, paste0(amu_updates_dir,'/unused_data_for_amu.xlsx'))
writexl::write_xlsx(unused_records_without_date, paste0(amu_updates_dir,'/unused_data_missing_dates.xlsx'))
#after editing DDDs, use this; amu_dataset_comb2

cat(paste("Done... uncleaned data stored", (length(unique(unused_records_amu$uid))), 'records\n'))
message(paste("Done... uncleaned data stored", (length(unique(unused_records_amu$uid))), 'records'))


#prepping the dataset
#creating the analysis dataset
amu_prev <- rbind(amu_dataset1,
                  amu_dataset_inhibitors %>% select(names(amu_dataset1)),
                  amu_dataset_comb1%>% select(names(amu_dataset1)),
                  amu_dataset_comb2%>% select(names(amu_dataset1)),
                  amu_dataset_comb_r_valid%>% select(names(amu_dataset1))) %>%
  filter(!is.na(date_of_data_collection)) %>%
  left_join(ddd_ref_update %>% dplyr::select(antibiotic_molecules, Class, aware_cats) %>% distinct(antibiotic_molecules, .keep_all = T), by='antibiotic_molecules')


age_factors <- c("0-27 days", "28-364 days", "1-4 years","5-9 years","10-14 years", "15 – 19 years", "20-24 years","25-59 years","60-99 years","100+ years")

amu_meta <- amu_raw_meta %>% select(facility,
                                    age,
                                    sex,
                                    ward_name,
                                    total_patients_in_ward,
                                    number_of_patients_included,
                                    uid, date_of_data_collection,
                                    patient_id
) %>%
  mutate(age_group = case_when(
    age >= 0       & age < 0.0767  ~ "0-27 days",        # ~27 days = 27/365
    age >= 0.0767  & age < 1       ~ "28-364 days",
    age >= 1       & age < 5       ~ "1-4 years",
    age >= 5       & age < 10      ~ "5-9 years",
    age >= 10      & age < 15      ~ "10-14 years",
    age >= 15      & age < 20      ~ "15–19 years",
    age >= 20      & age < 25      ~ "20-24 years",
    age >= 25      & age < 60      ~ "25-59 years",
    age >= 60      & age <= 99     ~ "60-99 years",
    age > 99   ~ "100+ years",
    TRUE ~ NA_character_
  ),
  age_group= factor(age_group, levels = age_factors)) %>%
  filter(!is.na(date_of_data_collection)) %>%
  dplyr::select(-date_of_data_collection) %>%
  filter(uid %in% amu_prev$uid)  #subsetting with confirmed antibiotic products


amu_prev_meta <- amu_prev %>%
  left_join(amu_meta, by='uid') %>%
  mutate(n=1,
         year=as.numeric(format(date_of_data_collection,'%Y')),
         month=as.numeric(format(date_of_data_collection,'%m')),
         y_month=format(as.Date(date_of_data_collection), "%Y-%m"),
         y_quarters=ifelse(month<=3, paste(year,'Q1'),
                           ifelse(month<=6 & month>3, paste(year,'Q1'),
                                  ifelse(month<=9 & month>6, paste(year,'Q1'),
                                         ifelse(month<=12 & month>9, paste(year,'Q1'),
                                                NA)))),
         y_month_date=ifelse(nrow(.)>0, as.Date(as.yearmon(y_month)),as.character(y_month)),
         aware_cats=ifelse(is.na(aware_cats), 'Uncategorized',aware_cats),
         Class=str_wrap(Class, width = 20)) %>%
  rename(antibiotic_class='Class',
         adm_route='route')


#Use frequency by other AMU variables
amu_vars_1 <- c("age_group" ,"antibiotic_class", "ward_name" , "sex", "indication_type", 'adm_route', "antibiotic_prescriber_type" ,
                "antibiotic_treatment_type", "antibiotic_guidelines_compliance", "antibiotic_oral_switch")


analysis_options <- amu_prev_meta %>% select(amu_vars_1) %>% summarise(across(everything(), ~ list(unique(.x)))) %>%
  pivot_longer(everything(),
               names_to = "variables_for_analysis",
               values_to = "options_in_dataset") %>%
  unnest(options_in_dataset) %>%
  mutate(user_standardized_options='')

cat(paste("Prepping done..."))
message(paste("Prepping done..."))

