cat('Beginning the analysis...\n')
cat(paste("Analysis proceeding with",(length(unique(amc_r1$uid))- length(unique(amc_dataset_comb2$uid))), 'records\n'))
message('Beginning the analysis...')
message(paste("Analysis proceeding with",(length(unique(amc_r1$uid))- length(unique(amc_dataset_comb2$uid))), 'records'))


#calculating the consumption (DiD)
amc1 <- amc_dataset1 %>% left_join(ddd_ref %>%
                                     dplyr::select(-uid,-antibiotic_names, -antibiotic_molecules, -antibiotic_names_orig_order) %>%
                                     distinct(name_route, .keep_all=T), by='name_route') %>%
  # filter(!is.na(`ATC code`)) %>%
  mutate(strength_unit=tolower(strength_unit),
         strength_unit=gsub('[-+_/,& )(;]','_',trimws(strength_unit)),
         strength_unit=str_split_i(strength_unit, '_',1),

         total=quantity*as.numeric(strength_val)*as.numeric(pack_size_val)/as.numeric(strength_div_factor_val),
         total_g=ifelse(strength_unit=='mg', total/1000,
                        ifelse(strength_unit=='mcg', total/1000000,
                        ifelse(strength_unit=='ml', total/1000,
                               ifelse(strength_unit=='g', total/1,
                                      NA)))),
         ddd_equivalent=total_g/as.numeric(DDD),
         year=as.numeric(format(date,'%Y')),
         y_month=format(as.Date(date), "%Y-%m"),
         y_month_date=as.Date(as.yearmon(y_month)),
         aware_cats=ifelse(is.na(aware_cats), 'Uncategorized',aware_cats)) #%>%

#continue here
#ddd for the inhibitors
amc2 <- amc_dataset_inhibitors %>% left_join(ddd_ref %>%
                                               dplyr::select(-uid,-antibiotic_names, -antibiotic_molecules, -antibiotic_names_orig_order) %>%
                                               distinct(name_route, .keep_all=T), by='name_route') %>%
  # filter(!is.na(`ATC code`)) %>%
  mutate(#strength_unit=str_split_i(strength_unit_r , "(/)", 1), #strength unit
    strength_unit=str_split_i(strength_num_factor,' ',1), #strength unit

         strength_div_factor1=str_split_i(strength_unit_r , "(/)",2),
         strength_div_factor_val1=str_split_i(strength_div_factor1 , "([A-Za-z]+)",1),
         strength_div_factor_val=ifelse(!is.na(strength_div_factor_val1), strength_div_factor_val1, strength_div_factor_val),
         strength_div_factor_val=ifelse(is.na(strength_div_factor_val)|strength_div_factor_val=='' , 1, strength_div_factor_val),
         strength_unit=tolower(strength_unit),

         strength_val1=gsub('[-+_/,& )(]','_',trimws(strength_val)),
         strength_val1=str_split_i(strength_val1, '_',1),
         total=quantity*as.numeric(strength_val1)*as.numeric(pack_size_val)/as.numeric(strength_div_factor_val),
         total_g=ifelse(strength_unit=='mg', total/1000,
                        ifelse(strength_unit=='mcg', total/1000000,
                               ifelse(strength_unit=='ml', total/1000,
                                      ifelse(strength_unit=='g', total/1,
                                             NA)))),
         ddd_equivalent=total_g/as.numeric(DDD),
         year=as.numeric(format(date,'%Y')),
         y_month=format(as.Date(date), "%Y-%m"),
         y_month_date=as.Date(as.yearmon(y_month)),
         aware_cats=ifelse(is.na(aware_cats), 'Uncategorized',aware_cats))


##the combination prods
amc3_1 <- amc_dataset_comb1 %>% filter(antibiotic_order==1) %>%
  left_join(ddd_ref %>%
              dplyr::select(-uid,-antibiotic_names, -antibiotic_molecules, -antibiotic_names_orig_order) %>%
              distinct(name_route, .keep_all=T), by='name_route') %>%
  mutate(#strength_unit=str_split_i(strength_unit_r , "(/)", 1), #strength unit
    strength_unit=str_split_i(strength_num_factor,' ',1), #strength unit
         strength_div_factor1=str_split_i(strength_unit_r , "(/)",2),
         strength_div_factor_val1=str_split_i(strength_div_factor1 , "([A-Za-z]+)",1),
         strength_div_factor_val=ifelse(!is.na(strength_div_factor_val1), strength_div_factor_val1, strength_div_factor_val),
         strength_div_factor_val=ifelse(is.na(strength_div_factor_val)|strength_div_factor_val=='' , 1, strength_div_factor_val),
         strength_unit=tolower(strength_unit),
         strength_unit=gsub('[-+_/,& )(;]','_',trimws(strength_unit)),
         strength_unit=str_split_i(strength_unit, '_',1),
    strength_unit=if_else(
      str_detect(strength_unit, "[0-9]+"),
      str_split_i(strength_unit, "([0-9]+)", 2),
      strength_unit
    ), #strength unit


         strength_val1=gsub('[-+_/,& )(;]','_',trimws(strength_val)),
         strength_val1=str_split_i(strength_val1, '_',1),
         total=quantity*as.numeric(strength_val1)*as.numeric(pack_size_val)/as.numeric(strength_div_factor_val),
         total_g=ifelse(strength_unit=='mg', total/1000,
                        ifelse(strength_unit=='mcg', total/1000000,
                               ifelse(strength_unit=='ml', total/1000,
                                      ifelse(strength_unit=='g', total/1,
                                             NA)))),
         ddd_equivalent=total_g/as.numeric(DDD),
         year=as.numeric(format(date,'%Y')),
         y_month=format(as.Date(date), "%Y-%m"),
         y_month_date=as.Date(as.yearmon(y_month)),
         aware_cats=ifelse(is.na(aware_cats), 'Uncategorized',aware_cats))


amc3_2 <- amc_dataset_comb1 %>% filter(antibiotic_order==2) %>%
  left_join(ddd_ref %>%
              dplyr::select(-uid,-antibiotic_names, -antibiotic_molecules, -antibiotic_names_orig_order) %>%
              distinct(name_route, .keep_all=T), by='name_route') %>%
  mutate(#strength_unit=str_split_i(strength_unit_r , "(/)", 1), #strength unit
         strength_unit=str_split_i(strength_num_factor,' ',1), #strength unit
         strength_div_factor1=str_split_i(strength_unit_r , "(/)",2),
         strength_div_factor_val1=str_split_i(strength_div_factor1 , "([A-Za-z]+)",1),
         strength_div_factor_val=ifelse(!is.na(strength_div_factor_val1), strength_div_factor_val1, strength_div_factor_val),
         strength_div_factor_val=ifelse(is.na(strength_div_factor_val)|strength_div_factor_val=='' , 1, strength_div_factor_val),
         strength_unit=tolower(strength_unit),
         strength_unit=gsub('[-+_/,& )(;]','_',trimws(strength_unit)),
         strength_unit=str_split_i(strength_unit, '_',1),

         strength_val1=gsub('[-+_/,&)(;]','_',trimws(strength_val)),
         strength_val1=str_split_i(strength_val1, '_',2),
         total=quantity*as.numeric(strength_val1)*as.numeric(pack_size_val)/as.numeric(strength_div_factor_val),
         total_g=ifelse(strength_unit=='mg', total/1000,
                        ifelse(strength_unit=='mcg', total/1000000,
                               ifelse(strength_unit=='ml', total/1000,
                                      ifelse(strength_unit=='g', total/1,
                                             NA)))),
         ddd_equivalent=total_g/as.numeric(DDD),
         year=as.numeric(format(date,'%Y')),
         y_month=format(as.Date(date), "%Y-%m"),
         y_month_date=as.Date(as.yearmon(y_month)),
         aware_cats=ifelse(is.na(aware_cats), 'Uncategorized',aware_cats))

##To continue adding cleanup options
#ddd for combinatio (separate the strengths)

#combine the amc dataset

amc <- Reduce(bind_rows, list(amc1 %>% dplyr::select(antibiotic_names, route,antibiotic_molecules, total_g, DDD, ddd_equivalent, year,y_month_date),
                              amc2%>% dplyr::select(antibiotic_names, route,antibiotic_molecules, total_g, DDD, ddd_equivalent, year,y_month_date),
                              amc3_1%>% dplyr::select(antibiotic_names, route,antibiotic_molecules=antibiotic_names_x1, total_g, DDD, ddd_equivalent, year,y_month_date),
                              amc3_2%>% dplyr::select(antibiotic_names, route,antibiotic_molecules=antibiotic_names_x1, total_g, DDD, ddd_equivalent, year,y_month_date))) %>%
  filter(!is.na(ddd_equivalent))

#subset million units also

#calculate the did later after merging
#group_by(year, antibiotic_names) %>%
#summarise()
#mutate(did=ddd_equivalent*1000/365/pop,)
#filter(ddd_equivalent>0)


##by molecule

amc_dir_molecule <- file.path(cntry,"Results_AMC",'molecules')

if(!dir.exists(amc_dir_molecule)){dir.create(amc_dir_molecule, recursive = T)}


amc %>% group_by(year, antibiotic_molecules) %>%
  summarise(  tot_ddd=sum(ddd_equivalent)) %>%
  ungroup() %>%
  mutate(tot_did=tot_ddd*1000/365/pop) %>%
  arrange(desc(tot_ddd)) %>%
  #  mutate(did=tot_ddd*1000/365/pop) %>%
  group_by(year) %>%
  mutate(annual_totals=sum(tot_ddd),
         ddd_dist=round(tot_ddd*100/annual_totals,2),
         molecule=factor(antibiotic_molecules, levels=unique(.$antibiotic_molecules)),
         molecule=str_to_title(molecule)) %>%
  arrange(desc(ddd_dist)) -> molecule_temp


amc_cats_molecule <- Reduce( rbind, list(molecule_temp %>% filter(ddd_dist >= 1),
                                         molecule_temp %>% filter(ddd_dist < 1) %>%   ##creating Others if less than 1%
                                           summarise(ddd_dist=sum(ddd_dist),
                                                     tot_did=sum(tot_did),
                                                     tot_ddd=sum(tot_ddd)) %>%
                                           mutate(molecule='Others')) ) %>%
  arrange(desc(ddd_dist)) %>%
  mutate( molecule=factor(molecule, levels=unique(.$molecule)))




#for (y in unique(amc_cats_molecule$year)) {


#To add - by sector


plt_molecule_dist <- ggplot(amc_cats_molecule,# %>% filter(year==y),
                            aes(x=as.factor(year), y=ddd_dist, fill=molecule))+
  geom_bar(stat = 'identity', width = 0.5)+
  scale_fill_manual(values = my_colors) +
  labs(x='Year', y='DDD Distribution (%)')+
  geom_label_repel(aes(label=ddd_dist, group=molecule),
                   position = position_stack(vjust = 0.5),
                   min.segment.length = 0,
                   segment.color = NA,    # removes leader lines
                   box.padding = 0.2,     # reduce push
                   show.legend = FALSE)+
  theme_classic()+
  theme(legend.title = element_blank())

ggsave(paste0(amc_dir_molecule,'/','AMC_molecules_dist_DDD.png'),plt_molecule_dist, width=8, height=8, units="in", dpi=300)

#next the barplots for did

plt_molecule_tot <- ggplot(amc_cats_molecule#%>% filter(year==y)
                           , aes(x=as.factor(year), y=tot_did, fill=molecule))+
  geom_bar(stat = 'identity',  width=.8, position = position_dodge())+
  labs(x='', y='DDD/1000 Inhabitants/day')+
  scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(legend.title = element_blank())

ggsave(paste0(amc_dir_molecule,'/','AMC_molecules_DID.png'),plt_molecule_tot, width=8, height=8, units="in", dpi=300)



plt_molecule_sin <- ggplot(amc_cats_molecule#%>% filter(year==y)
                           , aes(x=molecule, y=tot_did, fill=as.factor(year)))+
  geom_bar(stat = 'identity',  width=.8, position = position_dodge())+
  labs(x='', y='DDD/1000 Inhabitants/day')+
  scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),legend.title = element_blank())

ggsave(paste0(amc_dir_molecule,'/','AMC_molecules_single_DID.png'),plt_molecule_sin, width=8, height=8, units="in", dpi=300)

#}

cat('Analysis by antibiotic molecules done...\n')
message('Analysis by antibiotic molecules done...')
#s_routees distribution

cat('Begining analysis of antibiotic consumption by route of administration...\n')
message('Begining analysis of antibiotic consumption by route of administration...')


amc_dir_route <- file.path(cntry,"Results_AMC",'routes')

if(!dir.exists(amc_dir_route)){dir.create(amc_dir_route, recursive = T)}



amc %>% group_by(year, route) %>%
  summarise(  tot_ddd=sum(ddd_equivalent)) %>%
  ungroup() %>%
  mutate(tot_did=tot_ddd*1000/365/pop) %>%
  arrange(desc(tot_ddd)) %>%
  #  mutate(did=tot_ddd*1000/365/pop) %>%
  group_by(year) %>%
  mutate(annual_totals=sum(tot_ddd),
         ddd_dist=round(tot_ddd*100/annual_totals,2),
         route=factor(route, levels=unique(.$route)),
         route=str_to_title(route)) %>%
  arrange(desc(ddd_dist)) -> s_route_temp


amc_cats_s_route <- Reduce( rbind, list(s_route_temp %>% filter(ddd_dist >= 1),
                                        s_route_temp %>% filter(ddd_dist < 1) %>%   ##creating Others if less than 1%
                                          summarise(ddd_dist=sum(ddd_dist),
                                                    tot_did=sum(tot_did),
                                                    tot_ddd=sum(tot_ddd)) %>%
                                          mutate(s_route='Others')) ) %>%
  arrange(desc(ddd_dist)) %>%
  mutate( s_route=factor(s_route, levels=unique(.$s_route)))




#for (y in unique(amc_cats_s_route$year)) {


plt_s_route_dist <- ggplot(amc_cats_s_route #%>% filter(year==y),
                           ,aes(x=as.factor(year), y=ddd_dist, fill=route))+
  geom_bar(stat = 'identity', width = 0.5)+
  scale_fill_manual(values = my_colors) +
  labs(x='Year', y='DDD Distribution (%)')+
  geom_label_repel(aes(label=ddd_dist, group=route),
                   position = position_stack(vjust = 0.5),
                   min.segment.length = 0,
                   segment.color = NA,    # removes leader lines
                   box.padding = 0.2,     # reduce push
                   show.legend = FALSE)+
  theme_classic()+
  theme(legend.title = element_blank())

ggsave(paste0(amc_dir_route,'/','AMC_s_routes_dist_DDD.png'),plt_s_route_dist, width=8, height=8, units="in", dpi=300)

#next the barplots for did

plt_s_route_tot <- ggplot(amc_cats_s_route#%>% filter(year==y)
                          , aes(x=as.factor(year), y=tot_did, fill=route))+
  geom_bar(stat = 'identity',  width=.8, position = position_dodge())+
  labs(x='', y='DDD/1000 Inhabitants/day')+
  scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(legend.title = element_blank())
ggsave(paste0(amc_dir_route,'/','AMC_s_routes_DID.png'),plt_s_route_tot, width=8, height=8, units="in", dpi=300)

#}

cat('Analysis of antibiotic consumption by route done...\n')
message('Analysis of antibiotic consumption by route done...')



##by totals

amc_dir_total <- file.path(cntry,"Results_AMC",'totals')

if(!dir.exists(amc_dir_total)){dir.create(amc_dir_total, recursive = T)}


amc %>% group_by(year) %>%
  summarise(  tot_ddd=sum(ddd_equivalent)) %>%
  ungroup() %>%
  mutate(tot_did=tot_ddd*1000/365/pop) %>%
  arrange(desc(tot_ddd)) -> totals_temp


pd <- position_dodge2(width = 0.9, preserve = "single")

plt_did_tot <- ggplot(totals_temp#%>% filter(year==y)
                      , aes(x=as.factor(year), y=tot_did), fill='dodgerblue')+
  geom_col(position = pd)+
  labs(x='', y='Total DDD/1000 Inhabitants/day')+
  #scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(legend.title = element_blank())

ggsave(paste0(amc_dir_total,'/','AMC_Total_DID.png'),plt_did_tot, width=8, height=8, units="in", dpi=300)


#}

cat('Total DiDs done...\n')
message('Total DiDs done...')

##Consumption trend
# amc_dir_trend <- file.path(cntry,"Results_AMC",'Trend')
#
# if(!dir.exists(amc_dir_trend)){dir.create(amc_dir_trend, recursive = T)}
#
#
# amc_trend <- amc %>% group_by(y_month_date, Class, year) %>%
#   summarise(  tot_ddd=sum(ddd_equivalent)) %>%
#   ungroup() %>%
#   mutate(tot_did=tot_ddd*1000/365/pop) %>%
#   arrange(desc(tot_ddd))


#Visuals

cat('Beginning analysis of trends...\n')
message('Beginning analysis of trends...')

#classes distribution
amc_dir_class <- file.path(cntry,"Results_AMC",'class')

if(!dir.exists(amc_dir_class)){dir.create(amc_dir_class, recursive = T)}

class_names <- read_excel('amc_resources/ab_class_updated.xlsx') %>%
  dplyr::select(antibiotic_names, Class, Category) %>% distinct(antibiotic_names, .keep_all = T)

rearrange_words <- function(x, sep = ",") {
  sapply(strsplit(x, sep), function(words) {
    words <- trimws(words)                  # remove leading/trailing spaces
    paste(sort(words), collapse = sep)
  })
}

amc$antibiotic_names <- rearrange_words(amc$antibiotic_names)

unclassified_abs <- amc %>% left_join(class_names, by=('antibiotic_names')) %>%
  filter(is.na(Class)) %>%
  dplyr::select(antibiotic_names, antibiotic_molecules,Class, Category) %>%
  distinct(antibiotic_names, .keep_all=T)


cat('Analysis of trends done...\n')
message('Analysis of trends done...')
