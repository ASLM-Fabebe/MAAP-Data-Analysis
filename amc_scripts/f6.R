#read in the user-validated file

cat('Beginning AMC class analysis...\n')
message('Beginning AMC class analysis...')

file_path <- file.path(amc_updates_dir,'updated_AMC_classes.xlsx')

if (file.exists(file_path)) {
  amc_class_updates <- read.xlsx(file_path)# or read.csv, fread, etc.
} else {
  amc_class_updates=unclassified_abs %>% filter(!is.na(Class))
  message("File not found — amc_class_updates")
}

amc_class_updates=bind_rows_match_classes(list(class_names,amc_class_updates%>%
                                                 mutate(updated_by='User'))) %>%
  filter(!is.na(Class)) %>%
  distinct(antibiotic_names, .keep_all = T)

#update the original_file
write_xlsx(amc_class_updates,'amc_resources/ab_class_updated_b.xlsx')


amc %>% left_join(amc_class_updates %>% mutate(antibiotic_names=tolower(antibiotic_names)) %>%
                    distinct(antibiotic_names, .keep_all=T),
                  by=('antibiotic_molecules'='antibiotic_names')) %>%
  filter(!is.na(Class)) %>%   ##continue here
  group_by(year, Class) %>%
  summarise(  tot_ddd=sum(ddd_equivalent)) %>%
  ungroup() %>%
  mutate(tot_did=tot_ddd*1000/365/pop,
         Class=str_wrap(Class, width = 20)) %>%
  arrange(desc(tot_ddd)) %>%
  #  mutate(did=tot_ddd*1000/365/pop) %>%
  group_by(year) %>%
  mutate(annual_totals=sum(tot_ddd),
         ddd_dist=round(tot_ddd*100/annual_totals,2),
        # class=factor(Class, levels=unique(.$Class))
        ) %>%
  arrange((ddd_dist)) %>%
  mutate( cum=cumsum(ddd_dist),
          midpoint=ddd_dist/2+lag(cum))-> class_temp


amc_cats_class <- bind_rows_match_classes(list(class_temp %>% filter(ddd_dist >= 1),
                                      class_temp %>% filter(ddd_dist < 1) %>%   ##creating Others if less than 1%
                                        summarise(ddd_dist=sum(ddd_dist),
                                                  tot_did=sum(tot_did),
                                                  tot_ddd=sum(tot_ddd)) %>%
                                        mutate(Class='Others',
                                              # class=as.factor(class),
                                               midpoint=0)) ) %>%
  arrange(desc(ddd_dist)) %>%
  mutate( class=factor(Class, levels=c(unique(.$Class)[unique(.$Class) !='Others'], 'Others'))) %>%
  filter(!is.na(class))

#AWARe
amc_cats_aware <- amc %>% left_join(amc_class_updates, by=('antibiotic_names')) %>%
  filter(!is.na(Class)) %>%   ##continue here
  group_by(year, Category) %>%
  summarise(  tot_ddd=sum(ddd_equivalent)) %>%
  ungroup() %>%
  mutate(tot_did=tot_ddd*1000/365/pop) %>%
  arrange(desc(tot_ddd)) %>%
  #  mutate(did=tot_ddd*1000/365/pop) %>%
  group_by(year) %>%
  mutate(annual_totals=sum(tot_ddd),
         dist=round(tot_ddd*100/annual_totals,2),
         aware_cats=factor(Category, levels=c(unique(.$Category)[unique(.$Category) !='Others'], 'Others'))) %>%
  arrange(desc(dist)) %>%
  filter(!is.na(aware_cats))


#for (y in unique(amc_cats_class$year)) {



plt_class_dist <- ggplot(amc_cats_class #%>% filter(year==y)
                         , aes(x=as.factor(year), y=ddd_dist, fill=class))+
  geom_bar(stat = 'identity', width = 0.5)+
  scale_fill_manual(values = my_colors) +
  labs(x='Year', y='DDD Distribution (%)')+
  geom_label_repel(aes(label=ddd_dist, group=class),
                   position = position_stack(vjust = 0.5),
                   min.segment.length = 0,
                   segment.color = NA,    # removes leader lines
                   box.padding = 0.2,     # reduce push
                   show.legend = FALSE)+
  theme_classic()+
  theme(legend.title = element_blank())


write.csv(class_temp, paste0(amc_dir_class,'/','AMC_classes_dist.csv'))
ggsave(paste0(amc_dir_class,'/','AMC_classes_dist.png'),plt_class_dist, width=8, height=8, units="in", dpi=300)

#next the barplots for did

plt_class_tot <- ggplot(amc_cats_class#%>% filter(year==y)
                        , aes(x=as.factor(year), y=tot_did, fill=class))+
  geom_bar(stat = 'identity',  width=.8, position = position_dodge())+
  labs(x='', y='DDD/1000 Inhabitants/day')+
  scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(legend.title = element_blank())

ggsave(paste0(amc_dir_class,'/','AMC_classes.png'),plt_class_tot, width=8, height=8, units="in", dpi=300)


plt_class_sin <- ggplot(amc_cats_class#%>% filter(year==y)
                           , aes(x=class, y=tot_did, fill=as.factor(year)))+
  geom_col(position= position_dodge2(width = 0.9, preserve = "single"))+
  labs(x='', y='DDD/1000 Inhabitants/day')+
  scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),legend.title = element_blank())

ggsave(paste0(amc_dir_class,'/','AMC_single_classes.png'),plt_class_sin, width=8, height=8, units="in", dpi=300)




plt_aware_dist <- ggplot(amc_cats_aware %>% drop_na(aware_cats)#%>% filter(year==y)
                         , aes(x=as.factor(year), y=dist, fill=aware_cats))+
  geom_bar(stat = 'identity', width = 0.5)+
  #scale_fill_manual(values = my_colors) +
  labs(x='Year', y='DDD Distribution (%)')+
  scale_fill_manual(values=c('Access'='seagreen',
                             'Watch'='chocolate1',
                             'Reserve'='red',
                             'Uncategorized'='dodgerblue1',
                             'Not recommended'='purple3'
                             ))+
  geom_label_repel(aes(label=dist, group=aware_cats),
                   position = position_stack(vjust = 0.5),
                   min.segment.length = 0,
                   segment.color = NA,    # removes leader lines
                   box.padding = 0.2,     # reduce push
                   show.legend = FALSE)+
  theme_classic()+
  theme(legend.title = element_blank())

write.csv(amc_cats_aware, paste0(amc_dir,'/','AMC_AWARe_dist.csv'))
ggsave(paste0(amc_dir,'/','AMC_aware.png'),plt_aware_dist, width=8, height=8, units="in", dpi=300)


# plt_trend <- ggplot(amc_trend #%>% filter(year==y)
#                     , aes(x=y_month_date, y=log10(tot_ddd), color=Class))+
#   geom_line()+
#   labs(x='Year', y='DDDs (Log10)')+
#   theme_classic()+
#   theme(legend.title = element_blank(), legend.position = 'bottom')
#
# ggsave(paste0(amc_dir_trend,'/',y,'AMC_trends.png'),plt_trend, width=12, height=8, units="in", dpi=300)

#}

cat('Analysis Completed...\n')
message('Analysis Successfully Completed...')

cat(paste0('View all outputs in ', amc_dir))
message(paste0('View all outputs in ', amc_dir))


