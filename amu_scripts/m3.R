

cat('Hang tight Beginning clean up...\n')
message('Hang tight Beginning clean up...')


#read in the user-validated file
file_path <- file.path(amu_updates_dir,'matching_unclear_antibiotic_entries.xlsx')

if (file.exists(file_path)) {
  lookup_df <- read.xlsx(file_path)  # or read.csv, fread, etc.
} else {
  cat("File not found — keeping existing lookup_df\n")
}

# Split into matched and unmatched based on antibiotic_name NA
final_matched_df <- lookup_df[!is.na(lookup_df$antibiotic_name), ] %>%
  mutate(Source='User',
         # antibiotic_name=tolower(antibiotic_name),
         antibiotic_name=gsub('[-+_/,&]',',',tolower(antibiotic_name)),
         antibiotic_name=gsub('and ',',',tolower(antibiotic_name)),
         Verdict=as.character(Verdict)) %>%
  filter(grepl('have|corr', Verdict, ignore.case=T))

atcs_cleaned_update <-bind_rows_match_classes(list(atcs_cleaned, final_matched_df)) %>%
  distinct(original_entry, .keep_all = T)

#updating the atcs cleaned file
writexl::write_xlsx(atcs_cleaned_update, 'amc_resources/ab_molecules.xlsx')


#unmatched/non-antibiotic  entries
final_unmatched_df <- lookup_df[is.na(lookup_df$Verdict), ] %>%
  filter(nchar(original_entry)>3)

##
combined_matched_df <- bind_rows_match_classes(list(final_matched_df, matches_df)) %>%
  filter(!(antibiotic_name %in% tolower(inhibitors))) %>%
  distinct(original_entry, .keep_all = T)  #solving duplicates in the original_entry


#******************Mapping the values back
# Make a named vector for mapping: names are original entries, values are antibiotic names
map_vec <- setNames(combined_matched_df$antibiotic_name, combined_matched_df$original_entry)


# Function to map and replace values based on the lookup, preserving NAs and unmatched

map_antibiotics <- function(vec, mapping) {

  # For entries present in mapping, replace with mapped value; else keep original
  sapply(vec, function(x) {
    if (!is.na(x) && x %in% names(mapping)) {
      mapping[[x]]
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}


# Apply mapping on each namme* column
amu_r1_match_prep[n_cols] <- lapply(amu_r1_match_prep[n_cols], map_antibiotics, mapping = map_vec)

##creating the new ab_col retaining the inhibitors information
map_vec_all <- setNames(c(combined_matched_df$antibiotic_name, tolower(inhibitors)), c(combined_matched_df$original_entry, tolower(inhibitors)))

amu_r1_match_prep[n_cols] <- lapply(amu_r1_match_prep[n_cols], map_antibiotics, mapping = map_vec_all)



#************Returning the resolved dataframe
# Function to get antibiotic names from a vector and paste unique non-NA results

pattern_replace <- paste(other_meds_components, collapse = "|")


get_antibiotic_solved <- function(vec, lookup) {

  # Map entries to antibiotic names where possible
  mapped <- lookup[vec]

  # Remove NA and duplicates
  mapped_clean <- trimws(unique(mapped[!is.na(mapped)]))

  # Paste together or return NA if none
  if(length(mapped_clean) == 0) {
    return(NA_character_)
  } else {new_string = gsub(' ,', ',',paste(sort(mapped_clean), collapse = ","))   #order of combos in alphabetical order
  new_string = gsub(pattern_replace, "", new_string, ignore.case = TRUE)   #replacing the stabilizers and anasthetics
  return(
    paste(unique(trimws(unlist(strsplit(new_string, ",")))), collapse = ",")
  )  #for amu, have one without sort
  }
}


get_antibiotic_solved_original_order <- function(vec, lookup) {

  # Map entries to antibiotic names where possible
  mapped <- lookup[vec]
  # Remove NA and duplicates

  mapped_clean <- trimws(unique(mapped[!is.na(mapped)]))

  # Paste together or return NA if none
  if(length(mapped_clean) == 0) {
    return(NA_character_)
  } else {new_string = gsub(' ,', ',',paste((mapped_clean), collapse = ","))
  new_string = gsub(pattern_replace, "", new_string, ignore.case = TRUE)   #replacing the stabilizers and anasthetics
  return(
    paste(unique(trimws(unlist(strsplit(new_string, ",")))), collapse = ",")
  )  #for amu, have one without sort
  }
}


##for mapping onto the processed amu_r1_match_prep
map_vec_all_2 <- setNames(c(combined_matched_df$antibiotic_name, tolower(inhibitors)), c(combined_matched_df$antibiotic_name, tolower(inhibitors)))
map_vec_2 <- setNames(combined_matched_df$antibiotic_name, combined_matched_df$antibiotic_name)


# Apply per row  - we could introduce these 2 rows upfront in the dataset
amu_r1_match_prep$antibiotic_names <- apply(amu_r1_match_prep[n_cols], 1, get_antibiotic_solved, lookup = map_vec_all_2)
amu_r1_match_prep$antibiotic_molecules <- apply(amu_r1_match_prep[n_cols], 1, get_antibiotic_solved, lookup = map_vec_2)  #see tazobactam for example

#original combined antibiotics order (in the dataset)
amu_r1_match_prep$antibiotic_names_orig_order <- apply(amu_r1_match_prep[n_cols], 1, get_antibiotic_solved_original_order, lookup = map_vec_all_2)


##interrogating data with DDD information
amu_dataset_ed <-  amu_r1_match_prep%>% dplyr::select(-starts_with('x_name_part')) %>%
  mutate(name_route=paste0(antibiotic_names, '_',tolower(route)))

#records with defined ddd in the reference dataset
amu_dataset1 <- amu_dataset_ed%>%                   #CONTINUE HERE
  filter(name_route %in% ddd_ref$name_route)

#ddd of active comound is equal even in combination products (mostly in combinations with inhibitors)
amu_dataset_inhibitors <- subset(amu_dataset_ed, !(amu_dataset_ed$uid %in% amu_dataset1$uid)) %>%
  mutate(name_route=ifelse(!grepl(',', antibiotic_molecules) & grepl(',', antibiotic_names),
                           paste0(antibiotic_molecules, '_',tolower(route)),
                           name_route))%>%
  filter(name_route %in% ddd_ref$name_route)

#combination drugs
amu_dataset_comb <- subset(amu_dataset_ed, !(amu_dataset_ed$uid %in% c(amu_dataset1$uid, amu_dataset_inhibitors$uid))) %>%
  filter(!is.na(antibiotic_names))  #dropping all the NAs (products that could not be matched)


#breaking the combos
amu_dataset_comb1 <- amu_dataset_comb %>% mutate(antibiotic_names_x1=antibiotic_names_orig_order) %>%
  mutate(antibiotic_names_x1 = strsplit(antibiotic_names_orig_order, ",")) %>%
  #separate_rows(antibiotic_names_x1, sep = ",") %>%
  unnest_longer(antibiotic_names_x1, indices_to = "order") %>%
  rename(antibiotic_order = order) %>%
  mutate(name_route=paste0(antibiotic_names_x1, '_',tolower(route))) %>%
  filter(name_route %in% ddd_ref$name_route)

#not available
updates_comb <- amu_dataset_comb %>% mutate(antibiotic_names_x1=antibiotic_names_orig_order) %>%
  separate_rows(antibiotic_names_x1, sep = ",") %>%
  mutate(name_route=paste0(antibiotic_names_x1, '_',tolower(route))) %>%
  filter(!(name_route %in% ddd_ref$name_route)) %>%
  select(atc_level_name=antibiotic_names_x1, Adm.R =route) %>%
  distinct() %>%
  mutate(aware_cats= ' ',DDD=' ', Unit=' ') %>%
  filter(!is.na(atc_level_name))


#unused dataset for any ddd calculation
amu_dataset_comb_r <- subset(amu_dataset_comb, !(amu_dataset_comb$uid %in% c(amu_dataset1$uid, amu_dataset_inhibitors$uid, amu_dataset_comb1$uid)))


#
amu_dataset_comb_r_valid <- amu_dataset_comb_r %>% filter(antibiotic_molecules %in% str_split_i(ddd_ref$name_route, '_', 1))
#consider showing the incomplete data fir the user to clean/update like for the case of admin route


ddd_updates <- amu_dataset_comb_r %>% select(atc_level_name=antibiotic_molecules, Adm.R =route) %>%
  distinct() %>%
  mutate(aware_cats= '',DDD='', Unit='') %>%
  filter(!is.na(atc_level_name))

ddd_updates <- bind_rows_match_classes(list(ddd_updates, updates_comb)) %>%
  filter(!(atc_level_name %in% exclude_extra))

##the units to look up
units_options <- c('milligram', 'gram', 'millions of international units', 'international units', 'micrograms','unit dose', '' )

units_updates <- rbind(amu_dataset1 %>% distinct(antibiotic_unit_dose_units),
                       amu_dataset_inhibitors %>% distinct(antibiotic_unit_dose_units),
                       amu_dataset_comb1 %>% distinct(antibiotic_unit_dose_units)) %>%
  distinct(antibiotic_unit_dose_units) %>%
  mutate(standardized_units="") %>%
  rename(strength_unit_in_dataset=antibiotic_unit_dose_units) %>%
  filter(!(strength_unit_in_dataset %in% units_ref$strength_unit_in_dataset))





cat('Done!...')
message('Done!...')
