##

cat("Initializing updates on the dataset\n")
message("Initializing updates on the dataset\n")

##DDD updates
#read in the user-validated file
file_path <- file.path(amc_updates_dir,'DDD_information_updates.xlsx')

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
file_path <- file.path(amc_updates_dir,'strength_units_updates.xlsx')

if (file.exists(file_path)) {
  unit_updates <- read.xlsx(file_path)  # or read.csv, fread, etc.
} else {
  unit_updates=unit_updates
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
amc_dataset_comb2 <- amc_dataset_comb_r %>%  filter(name_route %in% ddd_ref$name_route)

unused_records <- subset(amc_r1, !(amc_r1$uid %in%  c(unique(amc_dataset1$uid),unique(amc_dataset_inhibitors$uid),
                                  unique(amc_dataset_comb1$uid), unique(amc_dataset_comb2$uid))))



writexl::write_xlsx(unused_records, paste0(amc_updates_dir,'/unused_amc_data.xlsx'))
#after editing DDDs, use this; amc_dataset_comb2

cat(paste("Done... uncleaned data stored", (length(unique(unused_records$uid))), 'records\n'))
message(paste("Done... uncleaned data stored", (length(unique(unused_records$uid))), 'records'))
