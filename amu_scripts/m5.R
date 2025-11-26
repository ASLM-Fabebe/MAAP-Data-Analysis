
#read in the user-validated file
file_path <- file.path(amu_updates_dir,'user_standardized_entries.xlsx')

if (file.exists(file_path)) {
  analysis_options <- read.xlsx(file_path)  # or read.csv, fread, etc.
} else {
  analysis_options=analysis_options
  cat("File not found — keeping existing lookup_df\n")
}

##integrate user defined options into the dataset
analysis_options <- analysis_options %>% filter(!is.na(user_standardized_options)) %>%
  mutate(id=paste0(row.names(.),variables_for_analysis))


for (i in analysis_options$id) {
  col_var=analysis_options$variables_for_analysis[analysis_options$id==i]
  col_opt=analysis_options$options_in_dataset[analysis_options$id==i]
  user_opt=analysis_options$user_standardized_options[analysis_options$id==i]

amu_prev_meta[[paste0(col_var)]][amu_prev_meta[[paste0(col_var)]]==paste0(col_opt)]=paste0(user_opt)

}
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
