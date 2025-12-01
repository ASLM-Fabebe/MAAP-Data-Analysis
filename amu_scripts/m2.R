
#updating the column names
cols_to_update <- read_excel(paste0(amu_updates_dir,"/select_amu_variables.xlsx")) %>%
  mutate(corresponding_variables=ifelse(is.na(corresponding_variables), 'not available', corresponding_variables))  #user canjust leave unavailable columns blank

#duplicated names - check if the user has provided duplicated names

cols_to_update_dup_check <- cols_to_update %>% filter(corresponding_variables!='not available')
dup_col <- cols_to_update_dup_check$corresponding_variables[duplicated((cols_to_update_dup_check$corresponding_variables))]

if (length(dup_col)>0) {
  message(paste0('The following columns are duplicated. Please correct before proceeding... ', dup_col))
}else{'no duplicates in your selected columns'}



#changing names
names(amu_raw)[names(amu_raw) %in% cols_to_update$corresponding_variables] <- cols_to_update$required_variables[match(names(amu_raw)[names(amu_raw) %in% cols_to_update$corresponding_variables], cols_to_update$corresponding_variables)]

#NAs for the unavailable cols
unavailable_cols=cols_to_update$required_variables[cols_to_update$corresponding_variables=='not available']

amu_raw[unavailable_cols]=NA

##prepping the metadata
amu_raw_meta <- amu_raw %>% dplyr::select(names(amu_raw)[!str_detect(names(amu_raw), '_opt')], uid) %>%
  mutate(age=ifelse(!is.na(age_yrs), as.numeric(age_yrs),
                    ifelse(!is.na(age_months), as.numeric(age_months)/12,
                           ifelse(!is.na(age_days), as.numeric(age_days)/365,
                                  NA))),
  date_of_data_collection=date_col_processing_vec(date_of_data_collection))%>%
    distinct()

ab_mix_hold <- list()

for (ab_mix in c('_opt1', '_opt2', '_opt3', '_opt4', '_opt5', '_opt6')) {

  if(length(names(amu_raw)[str_detect(names(amu_raw), paste0(ab_mix))])>0){

    x1<- amu_raw %>% dplyr::select(names(amu_raw)[str_detect(names(amu_raw), paste0(ab_mix))], uid, date_of_data_collection)   ##Note: loop over all the antibiotic options before doing the matches below

    #harmonize column names
    names(x1) <- gsub(paste0(ab_mix), '', names(x1))

    x2 <- x1 %>% mutate(antibiotic_names= str_split_i(antibiotic, "[0-9]",1),
                        antibiotic_copy=gsub('ampiclox','ampicillin and cloxacillin',tolower(antibiotic_names)),
                        antibiotic_copy=gsub('augmentin','amoxicillin and clavulanic acid',tolower(antibiotic_copy)),
                        antibiotic_copy=gsub('cotrimox','trimethoprim and sulfamethoxazole ',tolower(antibiotic_copy)),
                        antibiotic_copy=gsub('[-+_/,& )(;]','_s_',tolower(antibiotic_copy)),
                        antibiotic_copy=gsub(' and ','_s_',tolower(antibiotic_copy)),
                        antibiotic_unit_dose_units=tolower(antibiotic_unit_dose_units),
                        route=antibiotic_adm_route,
                        route=ifelse(tolower(antibiotic_adm_route)=='i','p',route))  #i for injection
    ab_mix_hold[[ab_mix]] <- x2

  }else{NULL}

}

amu_r1 <- do.call('rbind', ab_mix_hold) %>%
  filter(!is.na(antibiotic)) %>%
 # left_join(amu_raw_meta, by='uid') %>%    #rebuilding the dataset
  mutate(antibiotic_unit_dose_frequency= as.numeric(antibiotic_unit_dose_frequency),
         antibiotic_unit_dose_frequency= ifelse(antibiotic_unit_dose_frequency>6&antibiotic_unit_dose_frequency<30,
                                                24/antibiotic_unit_dose_frequency, antibiotic_unit_dose_frequency),#handling the 24h/half days entries
         antibiotic_therapy_start_date=date_col_processing_vec(antibiotic_therapy_start_date),
         antibiotic_therapy_end_date=date_col_processing_vec(antibiotic_therapy_end_date),
         #if end date therapy is missing, assume it is still active
         #antibiotic_therapy_end_date=ifelse(is.na(antibiotic_therapy_end_date), as.Date(date_of_data_collection), antibiotic_therapy_end_date),
         #if dose frequency is missing, assume it was administered once
         antibiotic_unit_dose_frequency=ifelse(is.na(antibiotic_unit_dose_frequency), 1, antibiotic_unit_dose_frequency),
         #if missed dose frequency is missing assume none
         antibiotic_missed_dose_frequency=ifelse(is.na(antibiotic_missed_dose_frequency), 0, antibiotic_missed_dose_frequency),
         #duration
         therapy_days=ifelse(is.na(antibiotic_therapy_end_date), (as.Date(date_of_data_collection)-as.Date(antibiotic_therapy_start_date)),(as.Date(antibiotic_therapy_end_date)-as.Date(antibiotic_therapy_start_date))),
         therapy_days=therapy_days+1  ##accounting for the initial day so 1-1=1
  )


#Resolving antibiotic names

# Determining the max number of parts after splitting
max_parts <- max(str_count(amu_r1$antibiotic_copy, "_s_"), na.rm = T) + 1

# Generating dynamic column names
col_names <- paste0("x_name_part", seq_len(max_parts))


amu_r1_match_prep <- amu_r1 %>% separate(antibiotic_copy, into = col_names, sep = "_s_", fill = "right")
amu_r1_match_prep[col_names] <- lapply(amu_r1_match_prep[col_names], trimws)  #trimming the whitespace


# Selecting columns starting with 'name'
n_cols <- grep("^x_name_part", names(amu_r1_match_prep), value = TRUE)

# Function to split into matched and unmatched
split_matches <- function(vec, dict_clean) {
  vec_clean <- trimws(tolower(vec))
  matched <- vec[vec_clean %in% dict_clean]
  unmatched <- vec[!vec_clean %in% dict_clean]
  list(
    matched = matched[!is.na(matched) & matched != ""],
    unmatched = unmatched[!is.na(unmatched) & unmatched != ""]
  )
}

# Apply to all 'name' columns
results_list <- lapply(n_cols, function(col) {
  split_matches(amu_r1_match_prep[[col]], antibiotics_mol_dict)
})

# Name by column
names(results_list) <- n_cols

# Combine all matches & unmatched into two unique vectors
all_matches   <- unique(unlist(lapply(results_list, `[[`, "matched")))
all_unmatched <- unique(unlist(lapply(results_list, `[[`, "unmatched")))

# Remove items in exclusion lists (case-insensitive)
final_unmatched <- all_unmatched[!all_unmatched %in% exclude_extra]

#if the data is perfect
if (length(final_unmatched)==0){
  final_unmatched=c('ampicillin')
}else{final_unmatched}


# Get antibiotic names (or NA)
cat('Hold tight! looking up potential matches in the Antibiotics database...\n')
message('Hold tight! looking up potential matches in the Antibiotics database...')

antibiotic_names <- ab_name(final_unmatched)

# Create lookup data frame for the unsure antibiotics
lookup_df <- data.frame(
  original_entry = final_unmatched,
  antibiotic_name = tolower(antibiotic_names),  ##might change when we implement the dynamic changes from user input. Would need creation of a new column

  source='amr_package',
  stringsAsFactors = FALSE) %>%
  mutate(antibiotic_name=gsub('/',',',antibiotic_name),
         Verdict='')


# Create df for exact matches (name = original)
matches_df <- data.frame(
  original_entry = sort(all_matches),
  antibiotic_name = (subset(atcs_cleaned, atcs_cleaned$original_entry %in% all_matches) %>% arrange(original_entry))$antibiotic_name,
  source='matched',
  stringsAsFactors = FALSE
)


# Combine both
# lookup_df <- lookup_df %>%
#   mutate()

cat('Done!...')
message('Done!...')
