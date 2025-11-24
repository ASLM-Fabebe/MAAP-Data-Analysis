
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
    age >= 15      & age < 20      ~ "15 – 19 years",
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

#plotting bars
pd <- position_dodge2(width = 0.9, preserve = "single")

#*Begin analysis
#*
#create the results directory
amu_dir <- file.path(cntry, "Results_AMU")

if(!dir.exists(amu_dir)){dir.create(amu_dir, recursive = T)}

#*Overall AMU prevalence % per facility
facility_prev_amu <- amu_prev_meta %>%
  arrange(desc(total_patients_in_ward)) %>%
  distinct(facility, ward_name, total_patients_in_ward, number_of_patients_included, year) %>% #picking the first reported number per month. Its a closer indicator to time overflow
  drop_na() %>%
  mutate(prev=number_of_patients_included/total_patients_in_ward) %>%
  filter(prev<=1) %>%   ##filtering out erroneous entries
  group_by(facility, year) %>%
  summarise(prev=round(mean(prev)*100,1))


if (nrow(facility_prev_amu)>0) {
overall_fac_plot <- ggplot(facility_prev_amu, aes(x = reorder(facility, -prev), y = prev, group = as.factor(year) )) +
  geom_col(width=.9, aes(fill = as.factor(year)))+  # Moved fill outside aes()
  labs(
    title = "Percentage of AMU across all facilities",
    x = "",
    y = "% of Patients on antibiotics"
  ) +
  scale_fill_manual(values = my_colors) +
    theme_classic() +
  theme(legend.position = 'right',legend.title = element_blank(),
        axis.text.x = element_text(angle = 0)) +
  coord_flip()
}else{p1=ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
  geom_point() +
  labs(x = "", y = "") +
  annotate("text",
           x = "a", y = 15,
           label = "Not enough data to generate this plot") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
}

ggsave(filename = paste0(amu_dir,"/1AMU_by_facility.png"), plot = overall_fac_plot, width = 6, height = 6, dpi = 300)


#Use frequency by other AMU variables
amu_vars_1 <- c("age_group" ,"antibiotic_class", "ward_name" , "sex", "indication_type", 'adm_route', "antibiotic_prescriber_type" ,
                "antibiotic_treatment_type", "antibiotic_guidelines_compliance", "antibiotic_oral_switch")

aware_colors <- c('Access'='seagreen',
                           'Watch'='chocolate1',
                           'Reserve'='red',
                           'Uncategorized'='dodgerblue1',
                           'Not recommended'='purple3'
)
#doing a loop over the variables
for (am_var in amu_vars_1) {

plot_facility_prev_amu(
  data = amu_prev_meta,
  x1 = am_var,
  my_colors = my_colors
)
}


##AWare
facility_aware_prev <- amu_prev_meta %>%
  distinct(
    facility, aware_cats,
    date_of_data_collection, patient_id, antibiotic_names,
    .keep_all = TRUE
  ) %>%
  group_by(facility, year) %>%
  mutate(tot_facility = sum(n)) %>%
  ungroup() %>%
  group_by(facility, aware_cats, year) %>%
  summarise(
    tot_facility = mean(tot_facility),
    ward_count   = n(),
    .groups = "drop") %>%
  mutate(prev = round(ward_count * 100 / tot_facility, 1),
         aware_cats=factor(aware_cats, levels = c('Access','Watch','Reserve', 'Not recommended','Uncategorized')))

write.csv(facility_aware_prev, paste0(amu_dir,"/AMU_by_AWaRe.csv"))


if (nrow(facility_aware_prev)>0) {

  p1 <- ggplot(
    facility_aware_prev,
    aes(
      x = reorder(facility, -prev),
      y = prev,
      group = as.factor(year)
    )
  ) +
    geom_col(width = .9, aes(fill = aware_cats)) +
    labs(
      #title = paste("Percentage of AMU across", x1),
      x = "",
      y = "% of antibiotics used"
    ) +
    scale_fill_manual(values = aware_colors) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 0)
    ) +
    coord_flip()
}else{p1=ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
  geom_point() +
  labs(x = "", y = "") +
  annotate("text",
           x = "a", y = 15,
           label = "Not enough data to generate this plot") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
}

ggsave(filename = paste0(amu_dir,"/AMU_by_AWaRe.png"), plot = p1, width = 6, height = 6, dpi = 300)


# AMU % by age group per facility
# AMU % by gender per facility
# Avg. antimicrobials per patient per facility
# AMU % by ATC class per facility
# AMU % by antimicrobial molecule per facility
# AMU % by AWaRe category per facility
# AMU % by administration route per facility


cat('Analysis Completed...\n')
message('Analysis Successfully Completed...')

cat(paste0('View all outputs in ', amu_dir))
message(paste0('View all outputs in ', amu_dir))
