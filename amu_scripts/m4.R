

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
