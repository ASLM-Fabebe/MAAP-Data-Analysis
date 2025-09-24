
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
famr_long_sir <- get_sir_interpr(df=famr_long)
famr_long_con <- get_con_interp(df=famr_long)


# Convert breakpoints to SIR ----------------------------------------------
message('Converting breakpoints to SIR ....')
cat('Converting breakpoints to SIR ....\n')

message('Please wait for this step to complete....')
cat('Please wait for this step to complete....\n')



amr_con <- convert2sir_fun(famr_long_con)

amr_sir <- convert2sir_fun(famr_long_sir)

# combine results

sir_outcomes_df <- dplyr::bind_rows(amr_con, amr_sir) %>%

  dplyr::filter(intrinsic_res_status=='FALSE') %>%   #drop the intrinsically resistant bug-drugs to not skew results

  dplyr::filter(interpreted_res!='NA')  #drop the UNINTERPRETABLE COMBOS FROM GUIDELINES


excluded_rec <- bind_rows(amr_con, amr_sir) %>%

  dplyr::filter(intrinsic_res_status=='TRUE'| interpreted_res=='NA')


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
loc_opts <- c('Inpatient','Outpatient','ICU or Critical care', 'Other')

loc_options=data.frame(my_dataset=unique(lkp_facility$`Patient Location Type`), options=c(rep('', length(unique(lkp_facility$`Patient Location Type`)))))


message('SIR interpretations complete ....')
cat('SIR interpretations complete  ....\n')










