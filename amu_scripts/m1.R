##load the packages
source('scripts/install_packages_packman.R')

#loading the functions
source('functions/amu_fun.R')

if (rstudioapi::isAvailable()) {
  folder_path <- rstudioapi::selectDirectory()
  print(folder_path)
} else {
  #  cat("Not running in RStudio.\n")
  folder_path=NA
}

#importing files both xls or xlsx
input_file <- list.files(
  folder_path,
  pattern = "^AMU.*\\.xls[x]?$",
  ignore.case = TRUE
)


amu_raw <- readxl::read_excel(file.path(folder_path,input_file)) %>%
    mutate(uid=1:nrow(.))


###matching the cols
choices1 <- c(names(amu_raw),'not available')  # Use your real variable



#prelim required columns
cols <- c("facility", "date_of_data_collection","patient_id",
          "ward_name","total_patients_in_ward","number_of_patients_included",
          "age_yrs", "age_months","age_days", 'sex','number_of_antibiotics',

          "indication_opt1", "indication_opt1_type","antibiotic_opt1","antibiotic_opt1_unit_dose",
          "antibiotic_opt1_unit_dose_from_combo","antibiotic_opt1_unit_dose_units",
          "antibiotic_opt1_unit_dose_frequency","antibiotic_opt1_adm_route",
          "antibiotic_opt1_therapy_start_date",  "antibiotic_opt1_therapy_end_date",
          "antibiotic_opt1_missed_dose_frequency", "antibiotic_opt1_sugical_proph_duration",
          "antibiotic_opt1_prescriber_type", "antibiotic_opt1_oral_switch",
          "antibiotic_opt1_guidelines_compliance","antibiotic_opt1_treatment_type",


          "indication_opt2", "indication_opt2_type","antibiotic_opt2","antibiotic_opt2_unit_dose",
          "antibiotic_opt2_unit_dose_from_combo","antibiotic_opt2_unit_dose_units",
          "antibiotic_opt2_unit_dose_frequency","antibiotic_opt2_adm_route",
          "antibiotic_opt2_therapy_start_date",  "antibiotic_opt2_therapy_end_date",
          "antibiotic_opt2_missed_dose_frequency", "antibiotic_opt2_sugical_proph_duration",
          "antibiotic_opt2_prescriber_type", "antibiotic_opt2_oral_switch",
          "antibiotic_opt2_guidelines_compliance","antibiotic_opt2_treatment_type",

          "indication_opt3", "indication_opt3_type","antibiotic_opt3","antibiotic_opt3_unit_dose",
          "antibiotic_opt3_unit_dose_from_combo","antibiotic_opt3_unit_dose_units",
          "antibiotic_opt3_unit_dose_frequency","antibiotic_opt3_adm_route",
          "antibiotic_opt3_therapy_start_date",  "antibiotic_opt3_therapy_end_date",
          "antibiotic_opt3_missed_dose_frequency", "antibiotic_opt3_sugical_proph_duration",
          "antibiotic_opt3_prescriber_type", "antibiotic_opt3_oral_switch",
          "antibiotic_opt3_guidelines_compliance","antibiotic_opt3_treatment_type",

          "indication_opt4", "indication_opt4_type","antibiotic_opt4","antibiotic_opt4_unit_dose",
          "antibiotic_opt4_unit_dose_from_combo","antibiotic_opt4_unit_dose_units",
          "antibiotic_opt4_unit_dose_frequency","antibiotic_opt4_adm_route",
          "antibiotic_opt4_therapy_start_date",  "antibiotic_opt4_therapy_end_date",
          "antibiotic_opt4_missed_dose_frequency", "antibiotic_opt4_sugical_proph_duration",
          "antibiotic_opt4_prescriber_type", "antibiotic_opt4_oral_switch",
          "antibiotic_opt4_guidelines_compliance","antibiotic_opt4_treatment_type",

          "indication_opt5", "indication_opt5_type","antibiotic_opt5","antibiotic_opt5_unit_dose",
          "antibiotic_opt5_unit_dose_from_combo","antibiotic_opt5_unit_dose_units",
          "antibiotic_opt5_unit_dose_frequency","antibiotic_opt5_adm_route",
          "antibiotic_opt5_therapy_start_date",  "antibiotic_opt5_therapy_end_date",
          "antibiotic_opt5_missed_dose_frequency", "antibiotic_opt5_sugical_proph_duration",
          "antibiotic_opt5_prescriber_type", "antibiotic_opt5_oral_switch",
          "antibiotic_opt5_guidelines_compliance","antibiotic_opt5_treatment_type",

          "indication_opt6", "indication_opt6_type","antibiotic_opt6","antibiotic_opt6_unit_dose",
          "antibiotic_opt6_unit_dose_from_combo","antibiotic_opt6_unit_dose_units",
          "antibiotic_opt6_unit_dose_frequency","antibiotic_opt6_adm_route",
          "antibiotic_opt6_therapy_start_date",  "antibiotic_opt6_therapy_end_date",
          "antibiotic_opt6_missed_dose_frequency", "antibiotic_opt6_sugical_proph_duration",
          "antibiotic_opt6_prescriber_type", "antibiotic_opt6_oral_switch",
          "antibiotic_opt6_guidelines_compliance","antibiotic_opt6_treatment_type")

empty_amu_df <- data.frame(required_variables=cols,
                           corresponding_variables=c(rep('',length(cols))))

#create a folder to hold the temporary files
amu_updates_dir <- file.path(folder_path, "analysis_updates")

#load_resources
source('amu_scripts/m1a.R')
