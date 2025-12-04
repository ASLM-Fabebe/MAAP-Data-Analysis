


# Functions ---------------------------------------------------------------

# Function 1: Convert breakpoints to SIR

convert2sir_fun <- function(df, n_cores = detectCores() - 1) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  id_vec <- df %>% pull(int_id)

  # Set up cluster
  cl1 <- makeCluster(n_cores)
  clusterEvalQ(cl1, {
    library(dplyr)
    library(AMR)
  })
  clusterExport(cl1, c("df"), envir = environment())

  res_list <- parLapply(cl1, id_vec, function(i) {

    df <- df %>% filter(!is.na(ab_name(drug_code)))  #drops the unrecognized antibiotic columns

    sub_df <- df %>% filter(int_id == i)

    id <- sub_df[["int_id"]]
    bacteria <- sub_df[["bacteria"]]
    drug_code <- as.ab(sub_df[["drug_code"]])
    test <- sub_df[["test_type"]]
    guideline <- sub_df[["guideline"]]
    value <- suppressWarnings(as.numeric(sub_df[["vals"]]))


    if (test == "mic") {
      if (is.na(value)) {
        value <- sub_df[["vals"]]
        interpreted <- as.character(as.sir(value, mo = bacteria, guideline = guideline, ab = drug_code))

      } else {
        interpreted <- as.character(as.sir(as.mic(value), mo = bacteria, guideline = guideline, ab = drug_code))


      }
    } else {
      if (is.na(value)) {
        value <- sub_df[["vals"]]
        interpreted <- as.sir(value, mo = bacteria, guideline = guideline, ab = drug_code)

      } else {
        interpreted <- as.sir(as.disk(value), mo = bacteria, guideline = guideline, ab = drug_code)



      }
    }

    #as_m1=as.character(interpreted)
    #as_d1=as.character(interpreted)

    intrinsic_status <- as.character(mo_is_intrinsic_resistant(bacteria, ab = drug_code))


    sub_df[["interpreted_res"]] <- as.character(interpreted)
    sub_df[["intrinsic_res_status"]] <- intrinsic_status
    sub_df[["breakpoints_available"]] <- ifelse(is.na(suppressWarnings(as.numeric(sub_df[["vals"]]))), as.integer(any(!is.na(c(as.character(as.sir(as.mic(3), mo = bacteria, guideline = guideline, ab = drug_code)),as.character(as.sir(as.disk(3), mo = bacteria, guideline = guideline, ab = drug_code)))))),
                                                2)  #if there are no breakpoints it will return zero


    sub_df
  })

  stopCluster(cl1)

  dplyr::bind_rows(res_list)
}

#------------------------------------------------------------------------------------------------

##Resistance calculation functions
indiv_ab_resistance <- function(df,path,path_par,org_res_dir_trends_par, org_res_dir_trends,...){
  hold_df <- df %>%
    filter(mo_organism==org_name & yr==y) %>%
    group_by(mo_organism, get(par_var_name), ab) %>%
    summarise(n=n(),
              r=sum(R),
              total_R = r/n  ) %>%
    mutate(ab_name=ab_name(ab),
           var_name=`get(par_var_name)`)

  write.csv(hold_df, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

  #plots
  ##preseving the empty cols for plotting equal bars
  hold_df <-  hold_df %>% filter(n>29) %>%
    filter(!is.na(var_name)) %>%
    ungroup() %>%
    complete(ab_name, var_name, fill = list(total_R = 0))%>%
    mutate(.is_missing = total_R == 0)   # flagging the ghosts

  ggsave(file.path(org_res_dir_par,paste0(cntry,'_',org_name,par,'_individual_abs','.png')), indiv_ab_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)

  #Analysis by period
  #yr

  if (y==unique(data_yrs$yr)[1]){
    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      group_by(mo_organism, get(par_var_name), ab, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             var_name=`get(par_var_name)`)

    write.csv(hold_df, file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,                # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    ggsave(file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,par,'_individual_abs','.png')), indiv_ab_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)

    #Quarters
    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      group_by(mo_organism, get(par_var_name), ab, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             var_name=`get(par_var_name)`)

    quarter_dir_par <- file.path(dirname(org_res_dir_trends_par),'Quarters')
    if(!dir.exists(quarter_dir_par)){dir.create(quarter_dir_par, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,                # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    ggsave(file.path(quarter_dir_par,paste0(cntry,'_',org_name,par,'_individual_abs','.png')), indiv_ab_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)


  }else{NULL}



  #trend end
}


overall_ab_resistance <- function(df,path,path_par, org_res_dir_trends, ...){

  if (y==unique(data_yrs$yr)[1]){   #this should be done once per session so it will be on the first year
    hold_df <- df %>%
      filter(mo_organism==org_name& yr==y) %>%
      #filter(ab %in% sel_abs) %>%
      group_by( ab) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab))

    write.csv(hold_df, file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_abs.csv')))
    #plot
    ggsave(file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')), overall_resistance_plot(hold_df),
           width=8, height=8, units="in", dpi=300)

    #Analysis by period
    #yr


    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      #filter(ab %in% sel_abs) %>%
      group_by( ab, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             n_lab=paste0('n=',n))

    write.csv(hold_df, file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs.csv')))
    #plot
    ggsave(file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs','.png')), overall_resistance_plot_trends(hold_df),
           width=8, height=8, units="in", dpi=300)

    #Quarters
    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      #filter(ab %in% sel_abs) %>%
      group_by( ab, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             n_lab=paste0('n=',n))

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs.csv')))
    #plot
    ggsave(file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')), overall_resistance_plot_trends_quarter(hold_df),
           width=8, height=8, units="in", dpi=300)
    #end trend

  }else{
    NULL
  }
}


indiv_ab_resistance_sau <- function(df, path,path_par,org_res_dir_trends_par, org_res_dir_trends,...){
  hold_df_mrsa <- df %>%
    filter(mo_organism==org_name & yr==y) %>%
    filter(ab %in% mrsa_abs) %>%
    group_by(mo_organism, get(par_var_name)) %>%
    summarise(n=n(),
              r=sum(R),
              total_R = r/n  ) %>%
    mutate(ab_name='MRSA',
           var_name=`get(par_var_name)`)


  write.csv(hold_df_mrsa, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_mrsa',par,'.csv')))


  #other combos
  hold_df_a <- df %>%
    filter(mo_organism==org_name& yr==y) %>%
    filter(ab %in% sel_abs) %>%
    group_by(mo_organism, get(par_var_name), ab) %>%
    summarise(n=n(),
              r=sum(R),
              total_R = r/n  ) %>%
    mutate(ab_name=ab_name(ab),
           var_name=`get(par_var_name)`)

  write.csv(hold_df_a, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

  hold_df <- bind_rows(hold_df_a, hold_df_mrsa)


  if (y==unique(data_yrs$yr)[1]){
    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0)   # flagging the ghosts

    ggsave(file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_overall_abs','.png')), indiv_ab_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)


    #*********analysis by period
    #*yr
    #*

    hold_df_mrsa <- df %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% mrsa_abs) %>%
      group_by(mo_organism, get(par_var_name), yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name='MRSA',
             var_name=`get(par_var_name)`)


    write.csv(hold_df_mrsa, file.path(org_res_dir_trends_par,paste0('year_',cntry,'_',org_name,'_mrsa',par,'.csv')))


    #other combos
    hold_df_a <- df %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% sel_abs) %>%
      group_by(mo_organism, get(par_var_name), ab, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             var_name=`get(par_var_name)`)

    write.csv(hold_df_a, file.path(org_res_dir_par,paste0('year_',cntry,'_',org_name,'_',par,'.csv')))

    hold_df <- bind_rows(hold_df_a, hold_df_mrsa)

    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    ggsave(file.path(org_res_dir_trends,paste0('year_',cntry,'_',org_name,'_overall_abs','.png')),
           indiv_ab_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)

    #Quarters
    hold_df_mrsa <- df %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% mrsa_abs) %>%
      group_by(mo_organism, get(par_var_name), yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name='MRSA',
             var_name=`get(par_var_name)`)

    quarter_dir_par <- file.path(dirname(org_res_dir_trends_par),'Quarters')
    if(!dir.exists(quarter_dir_par)){dir.create(quarter_dir_par, recursive = T)}


    write.csv(hold_df_mrsa, file.path(quarter_dir_par,paste0('year_',cntry,'_',org_name,'_mrsa',par,'.csv')))


    #other combos
    hold_df_a <- df %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% sel_abs) %>%
      group_by(mo_organism, get(par_var_name), ab, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             var_name=`get(par_var_name)`)

    quarter_dir_par <- file.path(dirname(org_res_dir_trends_par),'Quarters')
    if(!dir.exists(quarter_dir_par)){dir.create(quarter_dir_par, recursive = T)}

    write.csv(hold_df_a, file.path(quarter_dir_par,paste0('year_',cntry,'_',org_name,'_',par,'.csv')))

    hold_df <- bind_rows(hold_df_a, hold_df_mrsa)

    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}

    ggsave(file.path(quarter_dir,paste0('year_',cntry,'_',org_name,'_overall_abs','.png')),
           indiv_ab_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)

  }else{
    NULL
  }

  #Trend end
}

overall_ab_resistance_sau <- function(df, path,path_par,org_res_dir_trends,...){

  if (y==unique(data_yrs$yr)[1]){

    hold_df_mrsa <- df %>%
      filter(mo_organism==org_name& yr==y) %>%
      filter(ab %in% mrsa_abs) %>%
      group_by( mo_organism) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name='MRSA')


    write.csv(hold_df_mrsa, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'mrsa_overall_abs.csv')))

    #other combos
    hold_df_a <- df %>%
      filter(mo_organism==org_name & yr==y) %>%
      group_by( ab) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab))


    write.csv(hold_df_a, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_overall_abs.csv')))

    hold_df <- bind_rows(hold_df_a, hold_df_mrsa)

    #plots
    ggsave(file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_overall_abs','.png')), overall_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)


    #*****Analysis by period
    #*yr


    hold_df_mrsa <- df %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% mrsa_abs) %>%
      group_by( mo_organism, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name='MRSA')


    write.csv(hold_df_mrsa, file.path(org_res_dir_trends,paste0('year_',cntry,'_',org_name,'mrsa_overall_abs.csv')))

    #other combos
    hold_df_a <- df %>%
      filter(mo_organism==org_name) %>%
      group_by( ab,yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab))


    write.csv(hold_df_a, file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs.csv')))

    hold_df <- bind_rows(hold_df_a, hold_df_mrsa) %>%
      mutate(n_lab=paste0('n=',n))

    #plots
    ggsave(file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs','.png')),
           overall_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)

    #Quarters
    hold_df_mrsa <- df %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% mrsa_abs) %>%
      group_by( mo_organism, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name='MRSA')

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}

    write.csv(hold_df_mrsa, file.path(quarter_dir,paste0('year_',cntry,'_',org_name,'mrsa_overall_abs.csv')))

    #other combos
    hold_df_a <- df %>%
      filter(mo_organism==org_name) %>%
      group_by( ab,yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab))


    write.csv(hold_df_a, file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs.csv')))

    hold_df <- bind_rows(hold_df_a, hold_df_mrsa)%>%
      mutate(n_lab=paste0('n=',n))

    #plots
    ggsave(file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')),
           overall_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)
    #trend end
  }else{NULL
  }
}

#pathogen groups
indiv_ab_resistance_genus <- function(df,path,path_par,org_res_dir_trends_par, org_res_dir_trends, ...){
  hold_df <- df %>%
    mutate(mo_organism=genus ) %>%
    filter(ab %in% sel_abs) %>%
    filter(mo_organism==org_name & yr==y) %>%
    group_by(mo_organism, get(par_var_name), ab) %>%
    summarise(n=n(),
              r=sum(R),
              total_R = r/n  ) %>%
    mutate(ab_name=ab_name(ab),
           var_name=`get(par_var_name)`)

  write.csv(hold_df, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

  if (y==unique(data_yrs$yr)[1]){
    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0)   # flagging the ghosts

    ggsave(file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')), indiv_ab_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)

    #********Analysis by period
    #*yr

    hold_df <- df %>%
      mutate(mo_organism=genus ) %>%
      filter(ab %in% sel_abs) %>%
      filter(mo_organism==org_name ) %>%
      group_by(mo_organism, get(par_var_name), ab, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             var_name=`get(par_var_name)`)

    write.csv(hold_df, file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    ggsave(file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs','.png')),
           indiv_ab_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)

    #Quarters
    hold_df <- df %>%
      mutate(mo_organism=genus ) %>%
      filter(ab %in% sel_abs) %>%
      filter(mo_organism==org_name ) %>%
      group_by(mo_organism, get(par_var_name), ab, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             var_name=`get(par_var_name)`)

    quarter_dir_par <- file.path(dirname(org_res_dir_trends_par),'Quarters')
    if(!dir.exists(quarter_dir_par)){dir.create(quarter_dir_par, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir_par,paste0(cntry,'_',org_name,'_',par,'.csv')))

    #plots
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}

    ggsave(file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')),
           indiv_ab_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)
  }else{NULL
  }

  # trend end
}


overall_ab_resistance_genus <- function(df,path,path_par, org_res_dir_trends,...){

  if (y==unique(data_yrs$yr)[1]){

    hold_df <- df %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% sel_abs) %>%
      group_by( ab) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab))

    write.csv(hold_df, file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_abs.csv')))

    #plots
    ggsave(file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')),overall_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)

    #********Analysis by period
    #*Yr
    #*

    hold_df <- df %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% sel_abs) %>%
      group_by( ab, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             n_lab=paste0('n=',n))

    write.csv(hold_df, file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs.csv')))

    #plots
    ggsave(file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_abs','.png')),overall_resistance_plot_trends(hold_df),
           width=8, height=8, units="in", dpi=300)

    #Quarters
    hold_df <- df %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      filter(ab %in% sel_abs) %>%
      group_by( ab, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_name(ab),
             n_lab=paste0('n=',n))

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs.csv')))

    #plots
    ggsave(file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_abs','.png')),overall_resistance_plot_trends_quarter(hold_df),
           width=8, height=8, units="in", dpi=300)
  }else{NULL
  }

}




# Antibiotic class resistance ---------------------------------------------



#Individual pathogens
antibiotic_classes_res_indiv <- function(df,path,path_par,org_res_dir_trends_par, org_res_dir_trends,...) {

  # abx_classes <- unique(ab_group(abs_ref)) #returns the antibiotic classes

  hold_df <- df %>% #filter(ab_class %in% abx_classes) %>%  #turning this off for now
    filter(mo_organism==org_name & yr==y) %>%
    arrange(desc(R)) %>%
    distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
    group_by(mo_organism, get(par_var_name), ab_class) %>%
    summarise(n=n(),
              r=sum(R),
              total_R = r/n  ) %>%
    mutate(ab_name=ab_class,
           var_name=`get(par_var_name)`)

  write.csv(hold_df, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.csv')))

  #plot
  ##preseving the empty cols for plotting equal bars
  hold_df <-  hold_df %>% filter(n>29) %>%
    filter(!is.na(var_name)) %>%
    ungroup() %>%
    complete(ab_name, var_name, fill = list(total_R = 0))%>%
    mutate(.is_missing = total_R == 0)   # flagging the ghosts

  ggsave(file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.png')), indiv_ab_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)


  if (y==unique(data_yrs$yr)[1]){
    #calculating overall resistance
    hold_df <- df %>%
      filter(mo_organism==org_name & yr==y) %>%
      #filter(ab_class %in% abx_classes) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, ab_class) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class)

    write.csv(hold_df, file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_ab_classes.csv')))

    #plot
    ggsave(file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_ab_classes','.png')), overall_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)


    #***Analysis by periods

    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, get(par_var_name), ab_class, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             var_name=`get(par_var_name)`)

    write.csv(hold_df, file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.csv')))

    #plot
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    ggsave(file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.png')), indiv_ab_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)


    #calculating overall resistance
    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      #filter(ab_class %in% abx_classes) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, ab_class, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             n_lab=paste0('n=',n))

    write.csv(hold_df, file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_ab_classes.csv')))

    #plot
    ggsave(file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_ab_classes','.png')),
           overall_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)

    ##bY qUARTERS
    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, get(par_var_name), ab_class, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             var_name=`get(par_var_name)`)

    quarter_dir_par <- file.path(dirname(org_res_dir_trends_par),'Quarters')
    if(!dir.exists(quarter_dir_par)){dir.create(quarter_dir_par, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.csv')))

    #plot
    ##preseving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n

    ggsave(file.path(quarter_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.png')), indiv_ab_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)


    #calculating overall resistance
    hold_df <- df %>%
      filter(mo_organism==org_name) %>%
      #filter(ab_class %in% abx_classes) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, ab_class, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             n_lab=paste0('n=',n))

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_ab_classes.csv')))

    #plot
    ggsave(file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_ab_classes','.png')),
           overall_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)
  }else{NULL
  }

  #trend end

}


#Grouped pathogens
antibiotic_classes_res_grp <- function(df,path,path_par, org_res_dir_trends_par, org_res_dir_trends,...) {

  abx_classes <- unique(ab_group(abs_ref)) #returns the antibiotic classes

  hold_df <- df %>% #filter(ab_class %in% abx_classes) %>%
    mutate(mo_organism=genus) %>%
    filter(mo_organism==org_name & yr==y) %>%
    arrange(desc(R)) %>%
    distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
    group_by(mo_organism, get(par_var_name), ab_class) %>%
    summarise(n=n(),
              r=sum(R),
              total_R = r/n  ) %>%
    mutate(ab_name=ab_class,
           var_name=`get(par_var_name)`)

  write.csv(hold_df, file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.csv')))

  #plot
  ##preserving the empty cols for plotting equal bars
  hold_df <-  hold_df %>% filter(n>29) %>%
    filter(!is.na(var_name)) %>%
    ungroup() %>%
    complete(ab_name, var_name, fill = list(total_R = 0))%>%
    mutate(.is_missing = total_R == 0,
           n_lab=paste0('n=',n))   # flagging the ghosts

  ggsave(file.path(org_res_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.png')), indiv_ab_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)


  if (y==unique(data_yrs$yr)[1]){

    #calculating overall resistance
    hold_df <- df %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name & yr==y) %>%
      #filter(ab_class %in% abx_classes) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, ab_class) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             n_lab=paste0('n=',n))

    write.csv(hold_df, file.path(org_res_dir,paste0(cntry,'_',org_name,'_overall_ab_classes.csv')))

    #plot
    ggsave(file.path(org_res_dir,paste0(cntry,'_',org_name,cntry,'_overall_ab_classes','.png')), overall_resistance_plot(hold_df), width=8, height=8, units="in", dpi=300)

    #analysis by preriod
    #YR



    abx_classes <- unique(ab_group(abs_ref)) #returns the antibiotic classes

    hold_df <- df %>% #filter(ab_class %in% abx_classes) %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, get(par_var_name), ab_class,yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             var_name=`get(par_var_name)`)

    write.csv(hold_df, file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.csv')))

    #plot
    ##preserving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n


    ggsave(file.path(org_res_dir_trends_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.png')), indiv_ab_resistance_plot_trends(hold_df),
           width=8, height=8, units="in", dpi=300)


    #calculating overall resistance
    hold_df <- df %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      #filter(ab_class %in% abx_classes) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, ab_class, yr) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             n_lab=paste0('n=',n))

    write.csv(hold_df, file.path(org_res_dir_trends,paste0(cntry,'_',org_name,'_overall_ab_classes.csv')))

    #plot
    ggsave(file.path(org_res_dir_trends,paste0(cntry,'_',org_name,cntry,'_overall_ab_classes','.png')), overall_resistance_plot_trends(hold_df), width=8, height=8, units="in", dpi=300)

    #Quarters
    abx_classes <- unique(ab_group(abs_ref)) #returns the antibiotic classes

    hold_df <- df %>% #filter(ab_class %in% abx_classes) %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, get(par_var_name), ab_class,yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             var_name=`get(par_var_name)`)

    quarter_dir_par <- file.path(dirname(org_res_dir_trends_par),'Quarters')
    if(!dir.exists(quarter_dir_par)){dir.create(quarter_dir_par, recursive = T)}

    write.csv(hold_df, file.path(quarter_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.csv')))

    #plot
    ##preserving the empty cols for plotting equal bars
    hold_df <-  hold_df %>% filter(n>29) %>%
      filter(!is.na(var_name)) %>%
      ungroup() %>%
      complete(ab_name, var_name, fill = list(total_R = 0))%>%
      mutate(.is_missing = total_R == 0,   # flagging the ghosts
             n_lab=ifelse(is.na(n), '', paste0('n=',n)))   #indicating the n


    ggsave(file.path(quarter_dir_par,paste0(cntry,'_',org_name,'_',par,'ab_classes.png')), indiv_ab_resistance_plot_trends_quarter(hold_df),
           width=8, height=8, units="in", dpi=300)


    #calculating overall resistance
    hold_df <- df %>%
      mutate(mo_organism=genus) %>%
      filter(mo_organism==org_name) %>%
      #filter(ab_class %in% abx_classes) %>%
      arrange(desc(R)) %>%
      distinct(uid,mo_organism,ab_class, .keep_all = TRUE) %>% #remove duplicates for the Ab classes, depending on how vast the dataset is, could be updated to country, lab, etc
      group_by(mo_organism, ab_class, yr_quarter) %>%
      summarise(n=n(),
                r=sum(R),
                total_R = r/n  ) %>%
      mutate(ab_name=ab_class,
             n_lab=paste0('n=',n))

    quarter_dir <- file.path(dirname(org_res_dir_trends),'Quarters')
    if(!dir.exists(quarter_dir)){dir.create(quarter_dir, recursive = T)}


    write.csv(hold_df, file.path(quarter_dir,paste0(cntry,'_',org_name,'_overall_ab_classes.csv')))

    #plot
    ggsave(file.path(quarter_dir,paste0(cntry,'_',org_name,cntry,'_overall_ab_classes','.png')), overall_resistance_plot_trends_quarter(hold_df), width=8, height=8, units="in", dpi=300)

  }else{
    NULL
  }

}


# Plotting functions ------------------------------------------------------

overall_resistance_plot <- function(df) {
  if(nrow(df[df$n>29,])>0){
    plt_hold <- ggplot(df  %>% filter(n>29), aes(x=reorder(ab_name,total_R), y=total_R*100, group=ab_name))+
      geom_col(position = 'dodge', fill='dodgerblue')+
      labs(x='Antibiotic', y='Percent resistance')+
      theme_bw()+
      theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
            axis.title.x = element_text(size=11),
            plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
            strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='right')
    plt_hold
  }else{
    ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
      geom_point() +
      labs(x = "", y = "") +
      annotate("text",
               x = "a", y = 15,
               label = "Not enough data to generate this plot") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}


indiv_ab_resistance_plot <- function(df) {
  if(nrow(df[df$n>29,])>0){
    plt_hold <- ggplot(df, aes(x=ab_name, y=total_R*100, group=var_name, fill=var_name))+
      geom_col( position = position_dodge(width = 0.9))+
      #geom_text(aes( y=total_R*100,label=round(total_R*100,1)), position=position_dodge(width = 0.9), vjust=-.5)+
      labs(x=par, y='Percent resistance')+
      scale_fill_brewer(palette='Dark2',direction = -1)+
      theme_bw()+
      theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
            axis.title.x = element_text(size=11),
            #panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
            strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='bottom')
    plt_hold
  }else{
    ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
      geom_point() +
      labs(x = "", y = "") +
      annotate("text",
               x = "a", y = 15,
               label = "Not enough data to generate this plot") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}

indiv_ab_resistance_plot_trends = function(df) {
  if(nrow(df[df$n>29,])>0){
    pd <- position_dodge2(width = 0.9, preserve = "single")  # robust dodge
    # levels of your discrete fill
    lvls <- levels(factor(df$yr))
    # take the darker half of "Blues" and interpolate to the needed length
    pal_strong <- colorRampPalette(brewer.pal(9, "Blues")[4:9])(length(lvls))


    ggplot(df  %>% filter(n>29), aes(x=ab_name, y=total_R*100, fill=yr, group = yr))+
      facet_wrap(~var_name, ncol = 1)+
      geom_col(position = pd)+
      scale_fill_manual(values = setNames(pal_strong, lvls), name = "Period") +
      labs(x='', y='Percent resistance')+
      ylim(0,100)+
      #geom_text(aes(label=n_lab, group = yr), position=pd, angle=90, hjust=-.05)+
      theme_bw()+
      theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
            axis.title.x = element_text(size=11),
            plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
            strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='right')
  }else{
    ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
      geom_point() +
      labs(x = "", y = "") +
      annotate("text",
               x = "a", y = 15,
               label = "Not enough data to generate this plot") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}



overall_resistance_plot_trends <- function(df) {
  if(nrow(df[df$n>29,])>0){

    pd <- position_dodge2(width = 0.9, preserve = "single")  # robust dodge
    # levels of your discrete fill
    lvls <- levels(factor(df$yr))
    # take the darker half of "Blues" and interpolate to the needed length
    pal_strong <- colorRampPalette(brewer.pal(9, "Blues")[4:9])(length(lvls))


    ggplot(df  %>% filter(n>29), aes(x=ab_name, y=total_R*100, fill=yr, group = yr))+
      geom_col(position = pd)+
      scale_fill_manual(values = setNames(pal_strong, lvls), name = "Period") +
      labs(x='', y='Percent resistance')+
      ylim(0,100)+
      #geom_text(aes(label=n_lab, group = yr), position=pd, angle=90, hjust=-.05)+
      theme_bw()+
      theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
            axis.title.x = element_text(size=11),
            plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
            strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='right')

  }else{
    ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
      geom_point() +
      labs(x = "", y = "") +
      annotate("text",
               x = "a", y = 15,
               label = "Not enough data to generate this plot") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}


#Quarters
indiv_ab_resistance_plot_trends_quarter = function(df) {
  if(nrow(df[df$n>29,])>0){

    pd <- position_dodge2(width = 0.9, preserve = "single")  # robust dodge
    # levels of your discrete fill
    lvls <- levels(factor(df$yr_quarter))
    # take the darker half of "Blues" and interpolate to the needed length
    pal_strong <- colorRampPalette(brewer.pal(9, "Blues")[4:9])(length(lvls))


    ggplot(df  %>% filter(n>29), aes(x=ab_name, y=total_R*100, fill=yr_quarter, group = yr_quarter))+
      facet_wrap(~var_name, ncol = 1)+
      geom_col(position = pd)+
      scale_fill_manual(values = setNames(pal_strong, lvls), name = "Period") +
      labs(x='', y='Percent resistance')+
      ylim(0,100)+
      #geom_text(aes(label=n_lab, group = yr_quarter), position=pd, angle=90, hjust=-.05)+
      theme_bw()+
      theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
            axis.title.x = element_text(size=11),
            plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
            strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='right')

  }else{
    ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
      geom_point() +
      labs(x = "", y = "") +
      annotate("text",
               x = "a", y = 15,
               label = "Not enough data to generate this plot") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}



overall_resistance_plot_trends_quarter <- function(df) {
  if(nrow(df[df$n>29,])>0){

    pd <- position_dodge2(width = 0.9, preserve = "single")  # robust dodge
    # levels of your discrete fill
    lvls <- levels(factor(df$yr_quarter))
    # take the darker half of "Blues" and interpolate to the needed length
    pal_strong <- colorRampPalette(brewer.pal(9, "Blues")[4:9])(length(lvls))


    ggplot(df  %>% filter(n>29), aes(x=ab_name, y=total_R*100, fill=yr_quarter, group = yr_quarter))+
      geom_col(position = pd)+
      scale_fill_manual(values = setNames(pal_strong, lvls), name = "Period") +
      labs(x='', y='Percent resistance')+
      ylim(0,100)+
      #geom_text(aes(label=n_lab, group = yr_quarter), position=pd, angle=90, hjust=-.05)+
      theme_bw()+
      theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
            axis.title.x = element_text(size=11),
            plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
            strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='right')

  }else{
    ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
      geom_point() +
      labs(x = "", y = "") +
      annotate("text",
               x = "a", y = 15,
               label = "Not enough data to generate this plot") +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}





# Analysis functions ------------------------------------------------------


amr_individual_pathogens <-  function(xdf, org_res_dir,org_res_dir_par,org_name, abs_ref, cntry, par, par_var_name,org_res_dir_trends_par, org_res_dir_trends, ...){
  sel_abs <- subset(abs_ref,  abs_ref %in% xdf$ab)

  #subgroups
  environment(indiv_ab_resistance) <- environment()
  indiv_ab_resistance(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends_par, org_res_dir_trends,...)

  #calculating overall resistance
  environment(overall_ab_resistance) <- environment()
  overall_ab_resistance(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends,...)

  #class resistance
  environment(antibiotic_classes_res_indiv) <- environment()
  antibiotic_classes_res_indiv(xdf,org_res_dir, org_res_dir_par,org_res_dir_trends_par, org_res_dir_trends,...)

  if (exists("con") && inherits(con, "DBIConnection")) try(DBI::dbDisconnect(con), silent = TRUE)
  invisible(NULL)
}


#MRSA
mrsa_analysis <-  function(xdf,org_res_dir,org_res_dir_par,org_name, abs_ref, cntry, par, par_var_name, org_res_dir_trends_par, org_res_dir_trends,...){
  mrsa_abs <- c('FOX', 'MET', 'OXA')
  sel_abs_mrsa <- subset(mrsa_abs,  mrsa_abs %in% names(an_df))
  sel_abs <- subset(abs_ref,  abs_ref %in% names(an_df))

  #MRSA
  #subgroups
  environment(indiv_ab_resistance_sau) <- environment()
  indiv_ab_resistance_sau(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends_par, org_res_dir_trends,...)

  #calculating overall resistance
  #MRSA
  environment(overall_ab_resistance_sau) <- environment()
  overall_ab_resistance_sau(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends,...)

  #class resistance
  environment(antibiotic_classes_res_indiv) <- environment()

  antibiotic_classes_res_indiv(
    xdf,
    org_res_dir,          # path
    org_res_dir_par,      # path_par
    org_res_dir_trends_par,
    org_res_dir_trends,
    ...
  )
}

##
amr_pathogen_groups <-  function(xdf,org_res_dir,org_res_dir_par, org_name, abs_ref, cntry, par, par_var_name, org_res_dir_trends_par, org_res_dir_trends,...){
  sel_abs <- subset(abs_ref,  abs_ref %in% xdf$ab)

  #subgroups
  environment(indiv_ab_resistance_genus) <- environment()
  indiv_ab_resistance_genus(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends_par, org_res_dir_trends,...)

  #calculating overall resistance
  environment(overall_ab_resistance_genus) <- environment()
  overall_ab_resistance_genus(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends,...)

  #class resistance
  environment(antibiotic_classes_res_grp) <- environment()
  antibiotic_classes_res_grp(xdf,org_res_dir, org_res_dir_par, org_res_dir_trends_par, org_res_dir_trends,...)

  if (exists("con") && inherits(con, "DBIConnection")) try(DBI::dbDisconnect(con), silent = TRUE)
  invisible(NULL)
}


specimen_distribution <- function(df){

  spec_dir <- file.path(paste0(cntry,"/Results_AMR/GLASS"),'specimen_distribution')

  if(!dir.exists(spec_dir)){dir.create(spec_dir, recursive = T)}

  spec_plt_data <- df %>% distinct(uid, mo_organism, specimen_type, yr, .keep_all = T) %>%
    group_by(mo_organism, specimen_type, yr) %>%
    summarise(n_samples=n()) %>%
    mutate(n_lab=paste0('n=', n_samples))


  pd <- position_dodge2(width = 0.9, preserve = "single")  # robust dodge


  spec_dist_plot <- ggplot(spec_plt_data, aes(x=reorder(mo_organism, n_samples), y=n_samples, fill=specimen_type, group = specimen_type))+
    facet_wrap(~yr, ncol = 1)+
    geom_col(position = pd)+
    scale_fill_brewer(palette='Dark2',direction = 1)+
    labs(x='', y='Number of samples')+
    geom_text(aes(label=n_lab, group = yr), position=pd, angle=90, hjust=-.05)+
    ylim(0, max(spec_plt_data$n_samples)+.1*(max(spec_plt_data$n_samples)))+
    theme_bw()+
    theme(axis.text.x = element_text(size=11,angle=90,hjust=1,vjust=0),axis.text.y = element_text(size=10),axis.title = element_text(size=12),
          axis.title.x = element_text(size=11),
          plot.title = element_text(hjust = 0.5, face='bold'),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5,"cm"),
          strip.text = element_text(size=11),legend.text=element_text(size=8),legend.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position='right')

  ggsave(file.path(spec_dir,paste0(cntry,'_specimen_distribution','.png')), spec_dist_plot, width=8, height=8, units="in", dpi=300)

  spec_table <- spec_plt_data %>% dplyr::select(-n_lab)%>% pivot_wider(names_from = specimen_type, values_from = n_samples)
  write.csv(spec_table, file.path(spec_dir,paste0(cntry,'_specimen_table','.csv')))

  glass_labs <- df %>% distinct(r_id,specimen_type, .keep_all = T) %>%
    ungroup() %>%
    group_by(`Laboratory Name`) %>%
    summarise(Valid_GLASS_Records=n())

  write.csv(glass_labs, file.path(spec_dir,paste0(cntry,'_records_by_Labs','.csv')))
}

#------------------------------------------------------------------------------------------------------


