##Look up unclear entries
cat('checking your dataset...\n')
message('checking your dataset...')


###call and rename stuff

#create the results directory
amc_dir <- file.path(cntry, "Results_AMC")

if(!dir.exists(amc_dir)){dir.create(amc_dir, recursive = T)}


cols_to_update <- read_excel(paste0(amc_updates_dir,"/select_amc_variables.xlsx")) %>%
  mutate(Corresponding_variables=ifelse(is.na(Corresponding_variables), 'not available', Corresponding_variables))  #user canjust leave unavailable columns blank


# Perform renaming
names(amc_raw)[names(amc_raw) %in% cols_to_update$Corresponding_variables] <- cols_to_update$Required_variables[match(names(amc_raw)[names(amc_raw) %in% cols_to_update$Corresponding_variables], cols_to_update$Corresponding_variables)]


#creating blank space holders for the unavailable cols
unavailable_cols=cols_to_update$Required_variables[cols_to_update$Corresponding_variables=='not available']

amc_raw[unavailable_cols]=NA

if('date' %in% unavailable_cols) {amr$`date`=user_added_date }else{amr$`date`=amr$`date`} #fixes the lack od date column


#adding the is.na
amc_raw <- amc_raw %>% mutate(pack_size=ifelse(all(is.na(pack_size)),1, pack_size)) #solve when the quantity is already calculated, packsize is 1 so that the quantities are not overestimated



antibiotic_classes_amc <- c(
  "Aminoglycosides",
  "Amphenicols",
  "Antibiotics - Intestinal antiinfectives",
  "Carbapenems",
  "Combinations of drugs for treatment of tuberculosis",
  "Combinations of sulfonamides and trimethoprim, incl. derivatives",
  "First-generation cephalosporins",
  "Fluoroquinolones",
  "Fourth-generation cephalosporins",
  "Glycopeptide antibacterials",
  "Hydrazides",
  "Imidazole derivatives",
  "Lincosamides",
  "Macrolides",
  "Nitrofuran-derivatives",
  "Nitroimidazole derivatives",
  "Other agents against amoebiasis and other protozoal diseases",
  "Other drugs for treatment of tuberculosis",
  "Other cephalosporins",
  "Oxazolidinones",
  "Penicillins",
  "Polymyxins",
  "Rifamycins",
  "Second-generation cephalosporins",
  "Sulfonamides",
  "Tetracyclines",
  "Third-generation cephalosporins",
  "Triazole and tetrazole derivatives",
  "Tetracyclines",
  "Penicillins with extended spectrum",
  "Beta-lactamase sensitive penicillins",
  "Beta-lactamase resistant penicillins",
  "Beta-lactamase inhibitors",
  "Penicillins, incl. beta-lactamase inhibitors",
  "Monobactams",
  "Other cephalosporins and penems",
  "Trimethoprim and derivatives",
  "Short-acting sulfonamides",
  "Intermediate-acting sulfonamides",
  "Long-acting sulfonamides",
  "Streptogramins",
  "Streptomycins",
  "Other aminoglycosides",
  "Other quinolones",
  "Combinations of antibacterials",
  "Steroid antibacterials",
  "Nitrofuran derivatives",
  "Other antibacterials"
)

# Exclusion lists
inhibitors <- c(
  "Clavulanic acid", "Sulbactam", "Tazobactam", "Avibactam",'clavulanic',
  "Vaborbactam", "Relebactam", "Zidebactam", "Nacubactam",
  "PAβN", "D13-9001", "EDTA", "cilastatin", "betamipron", 'inhibitor'
)

other_meds_components <- c('benzathine',  #this is a stabilizer
                           'procaine') #this is an anesthetic)

other_non_antibiotics <- c(
  "diloxanide", "praziquantel", "antimalarial", "ors",'antibiotic','antifungal','cephalosporins',
  "i.v", "tabs", "cap", "syrup", "acid", "other",'oral', 'pack','plus', 'tablet', 'tablets', 'capsule','capsules'
)

# Combine exclusions
exclude_extra <- tolower(c(inhibitors, other_non_antibiotics))

#who ref
ddd_ref <- read_excel('amc_resources/ab_molecules_amc.xlsx') %>%
  mutate(name_route=tolower(name_route)) %>%
  distinct(name_route, .keep_all = T)

#atc molecules
atcs_cleaned <- read_excel('amc_resources/ab_molecules.xlsx') %>%
  mutate(original_entry=tolower(trimws(original_entry))) %>%
  distinct(original_entry, .keep_all = T)

# Clean ab dictionary
antibiotics_mol_dict <- trimws(tolower(atcs_cleaned$original_entry))

#strength units reference
units_ref <- read_excel('amc_resources/strength_units_updates_b.xlsx')


#load the dataset

#create the results directory
amc_dir <- file.path(cntry, "Results_AMC")
if(!dir.exists(amc_dir)){dir.create(amc_dir, recursive = T)}


#importing the test dataset (focus is on non-combinational drugs)
amc_test <- amc_raw %>%
  #amc_raw_sub %>%
  mutate(product_string=gsub("[0-9]", "", product)) %>%
  tidyr::extract("product", c("product1", "strength1"), "(\\D*)(\\d.*)", remove = F) %>%   #separating product names from strength
  tidyr::extract("pack_size", c("pack_size_unit1", "pack_size1"), "(\\D*)(\\d.*)") %>%   #separating packsize names from unit
 # mutate(product = if_else(is.na(product1), product, product1),
         mutate(product = product_string,

         #administration routes
         reported_route=route,
         route=ifelse(str_detect(tolower(route), 'oral'),'o',
                      ifelse(str_detect( tolower(route), 'parenteral'),'p',
                             tolower(route))),

         #If liquid formulation parameters are provided, paste the columns
         strength=ifelse(!is.na(strength_liquid_drugs), strength_liquid_drugs, strength),
         strength_unit=ifelse(!is.na(strength_liquid_drugs), strength_liquid_drugs, strength_unit),

         #populating strength, packsize and units if not provided
         strength=ifelse(is.na(strength), strength1, strength),

         # name_route=paste0(trimws(product), '_', route),
         strength_val=str_split_i(strength , "([A-Za-z]+)",1),   #strength integer (account for e+ of the millions of international units later)

         strength_unit_r=substr(strength , nchar(strength_val)+1, 50), #strength unit
         strength_unit1=str_split_i(strength_unit_r , "(/)", 1), #strength unit

         strength_unit=tolower(ifelse(strength_unit1==''|is.na(strength_unit1), strength_unit, strength_unit1)),

         #numerator for liquids
         #strength_num_factor=str_split_i(strength_unit , "(/)",1),
         #strength_num_factor_val=str_split_i(strength_num_factor , "([A-Za-z]+)",1),
         # strength_num_factor_val=ifelse((strength_num_factor_val==""),1,strength_num_factor_val),
         # strength_val=as.numeric(strength_val)*as.numeric(strength_num_factor_val),

         # strength_unit=str_split_i((substr(strength_num_factor , nchar(strength_num_factor_val)+1, 50)),' ',1))#, #strength unit

         # #div factor
         strength_div_factor=str_split_i(strength_unit_r , "(/)",2),
         strength_div_factor_val=str_split_i(strength_div_factor , "([A-Za-z]+)",1),
         strength_div_factor_val=ifelse(is.na(strength_div_factor_val),1 , strength_div_factor_val),
         strength_div_factor_unit=substr(strength_div_factor , nchar(strength_div_factor_val)+1, 50),   #strength unit
         pack_size_val=str_split_i(pack_size1 , "([A-Za-z]+)",1),   #packsize integer)
         #
         pack_unit1=substr(pack_size_unit1 , nchar(pack_size_val)+1, 50),
         #
         pack_unit=ifelse(is.na(pack_size_unit), pack_unit1, pack_size_unit),
         id=1:nrow(.))



#amc_ref=who_atc_ref

amc_r1 <- amc_test %>%
  mutate(
    antibiotic_copy=gsub('ampiclox','ampicillin and cloxacillin',tolower(product)),
    antibiotic_copy=gsub('augmentin','amoxicillin and clavulanic acid',tolower(antibiotic_copy)),
    antibiotic_copy=gsub('cotrimox','trimethoprim and sulfamethoxazole ',tolower(antibiotic_copy)),
    antibiotic_copy=gsub('[-+_/,& )(;]','_s_',tolower(antibiotic_copy)),
    antibiotic_copy=gsub(' and ','_s_',tolower(antibiotic_copy)),
    uid=1:nrow(.)                                    #adding an identifer to sort out the mixed molecules
  )

##taking care of date

amc_r1 <- amc_r1 %>%
  mutate(
    # pick source
    date_new = date_col_processing_vec(date))


# Determining the max number of parts after splitting
max_parts <- max(str_count(amc_r1$antibiotic_copy, "_s_"), na.rm = T) + 1

# Generating dynamic column names
col_names <- paste0("x_name_part", seq_len(max_parts))


amc_r1_match_prep <- amc_r1 %>% separate(antibiotic_copy, into = col_names, sep = "_s_", fill = "right")

amc_r1_match_prep[col_names] <- lapply(amc_r1_match_prep[col_names], trimws)  #trimming the whitespace


# Select columns starting with 'name'
n_cols <- grep("^x_name_part", names(amc_r1_match_prep), value = TRUE)

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
  split_matches(amc_r1_match_prep[[col]], antibiotics_mol_dict)
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



##*******************************Lookup erroneus entrries
# Assuming unmatched_filtered vector is already defined:

# Get antibiotic names (or NA)
cat('Hold tight! looking up potential matches in the Antibiotics database...\n')
message('Hold tight! looking up potential matches in the Antibiotics database...')

antibiotic_names <- ab_name(final_unmatched)

# Create lookup data frame
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
