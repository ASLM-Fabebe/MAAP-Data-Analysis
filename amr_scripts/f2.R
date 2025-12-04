
message('Loading functions ....')
cat('Loading functions ....')
# Load functions ----------------------------------------------------------

source(file.path("functions","amr_analysis_functions_01.R"))
source(file.path("functions","amr_analysis_functions_02.R"))


# Create output folder: Results -------------------------------------------


res_dir <- file.path(paste0(cntry,"/Results_AMR"))

if (!dir.exists(res_dir)){

  dir.create(res_dir, recursive = T)

}

date_var <- as.Date(date(), format = "%a %b %d %H:%M:%S %Y")



amr <- get_inputs_and_prelim_cleanup()

# Begin analysis here -----------------------------------------------------
message('Beginning analysis ....')
cat('Beginning analysis here ....')


lkp_demographics <- get_demographics(df=amr)

lkp_facility <- get_facilities_data(df=amr) %>%
  dplyr::select(r_id, 'Laboratory Name', 'Patient Department', "Patient Location Type")

lkp_specimens <- get_specimen_info(df=amr)

lkp_ast_ref <- get_guideline_info(amr)

amr_res <- get_test_results(df=amr)

famr_long <- pivot_abx_results(df=amr_res)

# separate breakpoints and SIR interpretations
sir_intp <- get_sir_interpr(df=famr_long)
con_intp <- get_con_interp(df=famr_long)

if(nrow(sir_intp)>0){famr_long_sir <- sir_intp%>%
  filter(!is.na(ab_name(drug_code)))}else{famr_long_sir =sir_intp[0, ]}

if(nrow(con_intp)>0){famr_long_con <- con_intp%>%
  filter(!is.na(ab_name(drug_code)))}else{famr_long_con =con_intp[0, ]}

  #famr_long_con <- get_con_interp(df=famr_long)%>% filter(!is.na(ab_name(drug_code)))


# Convert breakpoints to SIR ----------------------------------------------
message('Converting breakpoints to SIR ....')
cat('Converting breakpoints to SIR ....\n')

message('Please wait for this step to complete....')
cat('Please wait for this step to complete....\n')



#doing conversions in batches

chunk <- function(x, size=3000) {
  if (size <= 0) stop("size must be > 0")
  split(x, ceiling(seq_along(x) / size))
}

##mics and diameters
chunks <- chunk(famr_long_con$int_id, 3000)

if (length(unlist(chunks))==0){
  amr_con=NULL
}else{

message('Converting MICs and Disk measurements...')
cat('Converting MICs and Disk measurements...\n')

chunk_hold <- list()

for (ch in 1:length(chunks)){

  chunk_hold[[ch]] <- convert2sir_fun(famr_long_con %>% filter(int_id %in% unlist(chunks[ch])))

  message('Chunk....', ch, ' Interpretation completed, ', length(chunks)-ch, ' to go...')
  cat('Chunk....', ch, ' Interpretation completed, ', length(chunks)-ch, ' to go...\n')
}

amr_con <- do.call('rbind', chunk_hold)
}



#interpretations (SIR)
chunks <- chunk(famr_long_sir$int_id, 3000)

if (length(unlist(chunks))==0){
  amr_sir=NULL
}else{

  message('Validating your SIR interpretations...')
  cat('Validating your SIR interpretations...\n')

chunk_hold <- list()

for (ch in 1:length(chunks)){

  chunk_hold[[ch]] <- convert2sir_fun(famr_long_sir %>% filter(int_id %in% unlist(chunks[ch])))

  message('Chunk....', ch, ' Validation completed, ', length(chunks)-ch, ' to go...')
  cat('Chunk....', ch, ' Validation completed, ', length(chunks)-ch, ' to go...\n')
}

amr_sir <- do.call('rbind', chunk_hold)
}



# combine results
tmp <- purrr::compact(list(amr_con, amr_sir)) |>
  purrr::reduce(dplyr::bind_rows, .init = NULL)

if (!is.null(tmp) && nrow(tmp) > 0) {

sir_outcomes_df <- tmp %>%

  dplyr::filter(intrinsic_res_status=='FALSE') %>%   #drop the intrinsically resistant bug-drugs to not skew results

  dplyr::filter(breakpoints_available==1) %>% ##select  only valid tests even if SIRs are provided by the user

  dplyr::filter(interpreted_res!='NA') %>%   #drop the UNINTERPRETABLE COMBOS FROM GUIDELINES

  arrange(interpreted_res) %>%

  distinct(r_id, uid, organism,ab, .keep_all = T)

excluded_rec <- bind_rows(amr_con, amr_sir) %>%

  dplyr::filter(intrinsic_res_status=='TRUE'| interpreted_res=='NA')
} else {
  amc2 <- NULL
}
# get organism full names

lkp_organisms <- AMR::microorganisms %>% dplyr::select(mo,fullname)  #puls the entire list


# pivot wide

sir_outcomes_df_wide <- sir_outcomes_df %>%

  dplyr::select(-c(drug_code,int_id, vals, intrinsic_res_status)) %>%

  pivot_wider(names_from = "ab",

              values_from = "interpreted_res") %>%

  left_join(lkp_organisms, by=join_by("bacteria"=="mo")) %>%

  dplyr::select(r_id,uid,specimen_type,mo_organism=fullname,

                gramstain,test_type,guideline, everything())



#set up location lookup
#Use frequency by other AMU variables
amr_vars_1 <- c('Sex', 'specimen_type','location_type')


analysis_options <- sir_outcomes_df_wide %>% select(amr_vars_1) %>% summarise(across(everything(), ~ list(unique(.x)))) %>%
  pivot_longer(everything(),
               names_to = "variables_for_analysis",
               values_to = "options_in_dataset") %>%
  unnest(options_in_dataset) %>%
  mutate(user_standardized_options='')


#loc_opts <- c('Inpatient','Outpatient','ICU or Critical care', 'Other')

#loc_options=data.frame(my_dataset=unique(lkp_facility$`Patient Location Type`), options=c(rep('', length(unique(lkp_facility$`Patient Location Type`)))))

#setup specimen lookup
#
# spec_ops = sort(c(   "Blood",
#                      "Genital",
#                      "Respiratory",
#                      "Soft tissue & body fluids",
#                      "Stool",
#                      "Urine",
#                      "Other",
#                      "Unknown",
#                      "Cerebrospinal fluid",
#                      "Sputum", "Joint fluid"
# ))
#
# spec_all_options=data.frame(my_dataset=unique(lkp_facility$`Patient Location Type`), options=c(rep('', length(unique(lkp_facility$`Patient Location Type`)))))


message('SIR interpretations complete ....')
cat('SIR interpretations complete  ....\n')










