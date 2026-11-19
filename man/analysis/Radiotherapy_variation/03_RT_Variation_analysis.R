#' _____________________________________________________________________________
#' _____________________________________________________________________________

#### Variation in provider propensity to deliver frontline curative radiotherapy for lymphoma across England

#' _____________________________________________________________________________
#' _____________________________________________________________________________

rm(list=ls())
library(broom.mixed)
library(tidyverse)
library(readstata13)
library(FunnelPlotR)
library(COUNT)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(sjPlot)
library(lme4)
library(svglite)

#' _____________________________________________________________________________
#### Load data

raw_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/Analysis/Data/Raw/2023Update/"

analysed_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/Analysis/Data/Analysis/"

package_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/PHELymphRpackage/man/analysis/"

results_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/PHELymphRpackage/man/analysis/Radiotherapy_variation/results/"

load(paste0(raw_data_path, "tumour_reg.RData"))
df_lsoa <- tumour_reg %>% dplyr::select(pseudo_patientid, lsoa01_code) %>%
        rename(lsoa_code = lsoa01_code)

load(paste0(raw_data_path, "routes_diag.RData"))


# Gender = 1 male
## Load all HL and DLBCL patients, file written in NHS Trust identification file
## 39979, 9789HL and 30190 DLBCL patients
df_patients <- read.csv(file=paste0(package_data_path, 'Radiotherapy_variation/', "df_RT_analysis", '.csv'))

# But only 39955 on which we have a guess of org of care
df_orgs <- read.csv(file=paste0(package_data_path, 'Radiotherapy_variation/', "Patient_NHS_Trust_of_care", '.csv'))


#' _____________________________________________________________________________
#' Further exclusion based on organisation of care and generate analysis set

## NHS Trust names and code
# Many of the org codes in the patient records are old ones from old trust which have since merged

# df_orgcodes_orginal gives 2022 codes and names
df_orgcodes_orginal <- read.csv(paste0(package_data_path, "org_codes_original.csv"))

# df_orgcodes gives the possible older codes. E.g. BARTS HEALTH NHS TRUST is R1H in 2022 but includes R1H RGC RNH RNJ
# Mapping done in previous study https://doi.org/10.1111/bjh.19907Digital
df_orgcodes <- read.csv(paste0(package_data_path, "org_codes.csv"))
df_orgcodes_orginal <- read.csv(paste0(package_data_path, "org_codes_original.csv"))

#' ________
#' Analysis dataframe - create variables
#' Filter out those not diagnosed at acute trusts (i.e. mental health trust/ambulance - doesn't make sense)
#' Filter out DCOs
df_analysis <- left_join(df_orgs, df_patients, by = c("pseudo_patientid.x" = "pseudo_patientid.x")) %>%
        left_join(df_orgcodes, by = c("Org_of_care" = "trust.code")) %>%
        left_join(df_orgcodes_orginal, by = c("trust.name" = "trust.name")) %>%
        filter(trust.type == "ACUTE") %>%
        mutate(RT_binary = case_when(RT_intent_1 == "1L Curative" ~ 1, TRUE ~ 0)) %>%
        left_join(routes_diag,  by = c("pseudo_tumourid")) %>%
        filter(final_route != "DCO") %>%
        mutate(final_route_full = final_route) %>%
        mutate(final_route = factor(case_when(final_route =="" ~ "Unknown",
                                       final_route =="Emergency presentation" ~ "Emergency presentation",
                                       .default = "TWW/referral/outpatient/inpatient"), levels =c("TWW/referral/outpatient/inpatient","Emergency presentation", "Unknown"))) %>%
        filter(!is.na(chrl_tot_78_06)) %>%
        mutate(chrl_tot_78_06 = case_when(chrl_tot_78_06 > 3 ~ "3+",
                                          .default = as.character(chrl_tot_78_06))) %>%
        left_join(df_lsoa %>% distinct() %>% group_by(pseudo_patientid) %>% slice(1), by =c("pseudo_patientid.x" = "pseudo_patientid")) %>%
        dplyr:::select(pseudo_patientid.x, RT_binary, RT_intent_1, Org_of_care, trust.name, trust.code, gender, diagcat7,
               datdiag, diag_year, stage ,age_diag, imd, ethnicity, chrl_tot_78_06, final_route_full,final_route, lsoa_code) %>%
        mutate(age_diag_cat = case_when(age_diag > 17 & age_diag < 30~ "18-29",
                                        age_diag > 29 & age_diag < 40~ "30-39",
                                        age_diag > 39 & age_diag < 50~ "40-49",
                                        age_diag > 49 & age_diag < 60~ "50-59",
                                        age_diag > 59 & age_diag < 70~ "60-69",
                                        age_diag > 69 & age_diag < 80~ "70-79",
                                        .default = "80+")) %>%
        mutate(stage = as.factor(stage),
               stage = case_when(stage %in% c("5", "8", "99") ~ "Missing",.default = stage),
               stage_model = stage,
               stage_model = case_when( stage_model == "3" ~ "3/4",stage_model == "4" ~ "3/4",.default = stage_model),
               stage_model = factor(stage_model, levels = c("1","2","3/4","Missing"))) %>%
        fastDummies::dummy_cols(select_columns = c("stage_model"), remove_first_dummy = TRUE) %>%
        mutate(age_diag = as.numeric(age_diag)) %>%
        fastDummies::dummy_cols(select_columns = c("gender"), remove_first_dummy = TRUE) %>%
        rename(gender_Female = gender_2) %>%
        fastDummies::dummy_cols(select_columns = c("imd"), remove_first_dummy = TRUE) %>%
        fastDummies::dummy_cols(select_columns = c("chrl_tot_78_06"), remove_first_dummy = TRUE) %>%
        fastDummies::dummy_cols(select_columns = c("final_route"), remove_first_dummy = TRUE)

#' ________
#' Cut patients further
#### Look at Trust codes that are not on the list - may be private or older codes
df_NAorgcodes <- df_analysis %>%
        filter(is.na(trust.name)) %>%
        mutate(Org_of_care = as.character(Org_of_care)) %>%
        mutate(Org_of_care = as.factor(Org_of_care))

v_NAordcodes <- data.frame(levels(df_NAorgcodes$Org_of_care))

# Re,pve private
df_analysis <- df_analysis %>%
        filter(!is.na(trust.name))

## Need to remove certain trusts (and therefore patients)
df_analysis <- df_analysis %>%
        filter(!grepl(c("ALDER HEY|LIVERPOOL WOMEN'S NHS|HEREFORDSHIRE AND WORCESTERSHIRE|MIDLANDS PARTNERSHIP NHS|THE ROBERT JONES AND AGNES|THE WALTON CENTRE NHS FOUNDATION TRUST|VELINDRE NHS TRUST|THE ROYAL ORTHOPAEDIC HOSPITAL|MOORFIELDS EYE HOSPITAL|LIVERPOOL HEART AND|QUEEN VICTORIA HOSPITAL|ROYAL PAPWORTH|ROYAL NATIONAL ORTHOPAEDIC|GREAT ORMOND STREET"), trust.name))

## Counts of patients by trust
res_NHSTrust_patient <- df_analysis %>%
        group_by(diagcat7, trust.name) %>%
        count() %>%
        ungroup() %>%
        mutate(n = as.numeric(n)) %>%
        group_by(diagcat7) %>%
        summarise(median = median(n),
                  min = min(n),
                  max = max(n))

## 39712 patients now, chl:9717, dlbcl:29995
# df_analysis %>%
#         group_by(diagcat7 ) %>%
#         count()


#' _____________________________________________________________________________
#### Geographic data
## Cancer alliance map, Change it's CRS to match NHS trust data, Reorder, Make english boundary map

map_England_Cancer_alliances <- sf::st_read(paste0(package_data_path, "Geospatial_datasets/Cancer_Alliances/Cancer_Alliances_(July_2019)_EN_BFC.shp"), quiet = TRUE)
map_England_Cancer_alliances <- sf::st_transform(map_England_Cancer_alliances, crs = 4326)
map_England_Cancer_alliances <- map_England_Cancer_alliances %>%
        dplyr:::select(cal19nm, everything())
map_England_boundary <- sf::st_union(map_England_Cancer_alliances)

## NHS Trust map (point coordinates) and Filter to only those that are within English boundary
load(file = paste0(package_data_path, "Geospatial_datasets/map_points_nhs_trusts22.rda"))
points_nhs_trusts22 <- sf::st_filter(x = points_nhs_trusts22, y = map_England_boundary)

## LSOA centroids, Set up correct coordinate system, Filter to only those that are within English boundary
map_lsoa_centroids <- sf::st_read(paste0(package_data_path, "Geospatial_datasets/LSOA_Centroid_shape_file/LSOA_PopCentroids_EW_2021_V3.shp"), quiet = TRUE)
map_lsoa_centroids <- sf::st_transform(map_lsoa_centroids, crs = 4326)
map_lsoa_centroids <- sf::st_filter(x = map_lsoa_centroids, y = map_England_boundary)

## Load in 2011 to 2021 LSOA mapping. patient LSOA is 2011
mapping_LSOA2011_to_LSOA2021 <- read.csv(paste0(package_data_path, "Geospatial_datasets/LSOA_(2011)_to_LSOA_(2021)_to_Local_Authority_District_(2022)_Lookup_for_England_and_Wales_(Version_2).csv")) %>%
        dplyr::select(LSOA11CD,LSOA21CD)

mapping_LSOA2001_to_LSOA2011 <- read.csv(paste0(package_data_path, "Geospatial_datasets/LSOA2001_to_LSOA2011_mapping.csv"),
                                         na.strings = "") %>%
        dplyr::select(LSOA01CD,LSOA11CD)

# NHS trusts which have radiotherapy centres - use the coordinate as the place where the radiotherapy unit is
RT_centres_points_nhs_trusts22 <- points_nhs_trusts22 %>%
        filter(nhs_trust22_code %in% c("RGT","RN5","RA7","RWA","RTE","RDE","RNN","RK9","R0D","RJ1",
                                       "RYJ","RWF","RR8","RWE","RWD","RWH","RH5","RTD","RL4","RM1",
                                       "RAP","RNS","RX1","RTH","RGN","RHU","RRK","RF4","RHW","REF",
                                       "RTG","RH8","RAL","RPY","RXN","RXW","RA2","RYR","RD1","RAJ",
                                       "R1H","RBV","REN","RTR","RA9","RRV","RJE","RHM","RYG","RHQ",
                                       "RWP"))

## Add to analysis dataset
df_analysis_final <- df_analysis %>%
        # Persons location - a few cases of multiple locations due to LSOA code changes. cut those individuals
        left_join(mapping_LSOA2001_to_LSOA2011, by = c("lsoa_code" = "LSOA01CD")) %>%
        left_join(mapping_LSOA2011_to_LSOA2021, by = c("LSOA11CD")) %>%
        left_join(map_lsoa_centroids, by = c("LSOA21CD")) %>%
        distinct(pseudo_patientid.x, .keep_all = TRUE) %>%
        # location of each RT centre
        expand_grid(RT_centres_points_nhs_trusts22 %>%
                            dplyr::select(nhs_trust22_code, geometry ) %>%
                            rename(NHSTrust_geometry = geometry)) %>%
        rename(`LSOA geometry` = 'geometry') %>%
        mutate(temp = gsub('[c()]', '', NHSTrust_geometry)) %>%
        separate(col = temp, into = c('Longitude_trust', 'Latitude_trust'), sep = '\\,') %>%
        mutate(temp = gsub('[c()]', '', `LSOA geometry`)) %>%
        separate(col = temp, into = c('Longitude_LSOA', 'Latitude_LSOA'), sep = '\\,') %>%
        mutate(across(Longitude_trust:Latitude_LSOA, as.numeric)) %>%
        mutate(distance_km = round(geosphere::distHaversine(cbind(Longitude_LSOA, Latitude_LSOA), cbind(Longitude_trust, Latitude_trust)) /1000, 2    )) %>%
        filter(!is.na(distance_km))%>%
        group_by(pseudo_patientid.x) %>%
        slice(which.min(distance_km)) %>%
       mutate(distance_km_cat = case_when(distance_km > 0 & distance_km <= 10 ~ "0-10",
                                   distance_km > 10 & distance_km <= 20~ "10-20",
                                   distance_km > 20 & distance_km <= 30~ "20-30",
                                   distance_km > 30 & distance_km <= 40~ "30-40",
                                   .default = "40+"))


#' _________________________________________________________________________________
### Results - 3.1

# Patient characteristic tables
res_descriptives_n <- df_analysis_final %>%
        group_by(diagcat7) %>%
        summarise(total = n())

res_descriptives_age <- df_analysis_final %>%
        group_by(diagcat7) %>%
        summarise(median = median(age_diag),
                  iqr_75 = quantile(age_diag, 3/4),
                  iqr_25 = quantile(age_diag, 1/4))

res_descriptives_age_cat <- df_analysis_final %>%
        group_by(diagcat7, age_diag_cat) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        group_by(diagcat7) %>%
        mutate(all = sum(total)) %>%
        ungroup() %>%
        mutate(percentage = round((total / all) *100 ,1))

res_descriptives_sex <- df_analysis_final %>%
        group_by(diagcat7, gender) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        group_by(diagcat7) %>%
        mutate(total_diagcat = sum(total),
               percentage_diagcat = round((total/total_diagcat) *100,2))

res_descriptives_stage <- df_analysis_final %>%
        group_by(diagcat7, stage) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        group_by(diagcat7) %>%
        mutate(total_diagcat = sum(total),
               percentage_diagcat = round((total/total_diagcat) *100,1))

res_descriptives_IMD <- df_analysis_final %>%
        group_by(diagcat7, imd) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        group_by(diagcat7) %>%
        mutate(total_diagcat = sum(total),
               percentage_diagcat = round((total/total_diagcat) *100,1))

res_descriptives_CCI <- df_analysis_final %>%
        group_by(diagcat7, chrl_tot_78_06) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        group_by(diagcat7) %>%
        mutate(total_diagcat = sum(total),
               percentage_diagcat = (total/total_diagcat) *100)

res_descriptives_route <- df_analysis_final %>%
        mutate(final_route_full = case_when( final_route_full == "" ~ "Unknown",
                                             .default = final_route_full)) %>%
         group_by(diagcat7, final_route_full) %>%
        summarise(total = n()) %>%
        ungroup() %>%
        group_by(diagcat7) %>%
        mutate(total_diagcat = sum(total),
               percentage_diagcat = round((total/total_diagcat) *100,1))

# Proportion of patients recieving radiotherapy over time
res_prop_RT_over_time <- df_analysis_final %>%
        group_by(diagcat7, stage , diag_year, RT_binary) %>%
        summarise(count = n() ) %>%
        mutate(total = sum(count)) %>%
        mutate(prop = round((count / sum(count)) * 100, 0) ) %>%
        filter(RT_binary == 1) %>%
        mutate(label = case_when(diagcat7 == "CHL" ~ paste0(count," of ",total, " HL patients"),
                                 diagcat7 == "DLBCL" ~ paste0(count," of ",total, " DLBCL patients")) ) %>%
        mutate(diag_year = as.factor(diag_year))

res_descriptives_distance <- df_analysis_final %>%
        group_by(diagcat7) %>%
        summarise(median = median(distance_km),
                  iqr_75 = quantile(distance_km, 3/4),
                  iqr_25 = quantile(distance_km, 1/4))

res_descriptives_distance_cat <- df_analysis_final %>%
        mutate(distance_km_cat = case_when(distance_km > 0 & distance_km <= 10 ~ "0-10",
                                           distance_km > 10 & distance_km <= 20~ "10-20",
                                           distance_km > 20 & distance_km <= 30~ "20-30",
                                           distance_km > 30 & distance_km <= 40~ "30-40",
                                           .default = "40+")) %>%
                                           group_by(diagcat7, distance_km_cat) %>%
                                                   summarise(total = n()) %>%
                                                   ungroup() %>%
                                                   group_by(diagcat7) %>%
                                                   mutate(total_diagcat = sum(total),
                                                          percentage_diagcat = round((total/total_diagcat) *100,2))


# 1000 x 600
plot_fig1 <- ggplot(res_prop_RT_over_time %>% mutate(diagcat7 = case_when(diagcat7 == "CHL" ~ "cHL", .default = diagcat7 ),
                                                     diagcat7 = factor(diagcat7, levels = c("cHL", "DLBCL"))),
                    aes(x=diag_year, y= prop, group = stage, colour = stage , label = label)) +
        geom_line()+
        geom_point() +
        facet_wrap( facets = "diagcat7", nrow = 2 ) +
        scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100),labels = scales::label_percent(scale = 1)) +
        theme_bw(base_size = 16) +
        ggsci::scale_colour_jco(labels = c("I","II","III","IV","Missing")) +
        ylab("Percentage of newly diagnosied patients\nreceiving frontline curative intent radiotherapy") +
        xlab("Diagnosis year") +
        labs(color='Stage at diagnosis') +
        theme(legend.position="bottom",
              strip.background = element_rect(fill="white"))

ggsave(paste0( results_path, "Figure_1_prop_variation.png"), plot = plot_fig1, width = 19, height =19, units = "cm", device='png', dpi=600)

# Proportion of patients recieving radiotherapy by stage
res_prop_RT_bystage <- df_analysis_final %>%
        group_by(diagcat7, stage_model, RT_binary) %>%
        summarise(count = n()) %>%
        mutate(total = sum(count)) %>%
        mutate(prop = round((count / sum(count)) * 100, 0) ) %>%
        filter(RT_binary == 1) %>%
        mutate(label = case_when(diagcat7 == "CHL" ~ paste0(count," of ",total, " HL patients"),
                                 diagcat7 == "DLBCL" ~ paste0(count," of ",total, " DLBCL patients")) )

#' _____________________________________________________________________________
##__________ Results - 3.2 Institutional comparison

## Logistic regression to calculate expected odds of radiotherapy using logistic regression

## For Hl and DLBCL seperate - then correlation

df_analysis_final_cHL <- df_analysis_final %>% filter(diagcat7 == "CHL") %>%
        fastDummies::dummy_cols(select_columns = c("age_diag_cat"), remove_first_dummy = TRUE) %>%
        fastDummies::dummy_cols(select_columns = c("distance_km_cat"), remove_first_dummy = TRUE)

expected_model_HL <- glm(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                                 imd_2 + imd_3  + imd_4  + imd_5 +
                                 chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                                 `final_route_Emergency presentation` + final_route_Unknown + `age_diag_cat_30-39` + `age_diag_cat_40-49` +
                                 `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                                 `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                                 `distance_km_cat_40+`, family=binomial(link='logit'), data = df_analysis_final_cHL)
res_logistic_HL_AME <- tidy(expected_model_HL,conf.int=TRUE, exponentiate=TRUE, effects="fixed")

df_analysis_final_DLBCL <- df_analysis_final %>% filter(diagcat7 == "DLBCL") %>%
        mutate(age_diag_cat = case_when(age_diag > 17 & age_diag < 50~ "18-49",
                                        age_diag > 49 & age_diag < 60~ "50-59",
                                        age_diag > 59 & age_diag < 70~ "60-69",
                                        age_diag > 69 & age_diag < 80~ "70-79",
                                        .default = "80+"
        )) %>%
        fastDummies::dummy_cols(select_columns = c("age_diag_cat"), remove_first_dummy = TRUE) %>%
        fastDummies::dummy_cols(select_columns = c("distance_km_cat"), remove_first_dummy = TRUE)


expected_model_DLBCL <- glm(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                                    imd_2 + imd_3  + imd_4  + imd_5 +
                                    chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                                    `final_route_Emergency presentation` + final_route_Unknown +
                                    `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                                    `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                                    `distance_km_cat_40+`  , family=binomial(link='logit'), data = df_analysis_final_DLBCL)
res_logistic_DLBCL_AME <- tidy(expected_model_DLBCL,conf.int=TRUE, exponentiate=TRUE, effects="fixed")


df_analysis_final_cHL$expected <- predict(expected_model_HL, type = "response")
df_analysis_final_DLBCL$expected <- predict(expected_model_DLBCL, type = "response")

## For the funnel plot analysis this package is used.
#devtools::install_github("oizin/funnelplot")
library(funnelplot)

# outcome ~ covariates | cluster ID
f1 <- funnelModel(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                          imd_2 + imd_3  + imd_4  + imd_5 +
                          chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                          `final_route_Emergency presentation` + final_route_Unknown + `age_diag_cat_30-39` + `age_diag_cat_40-49` +
                          `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                          `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                          `distance_km_cat_40+`|
                          Org_of_care ,
                  pointTarget(limits = c(0.05,0.003), standardised = FALSE, normalApprox = TRUE), data = df_analysis_final_cHL)

# Note warnings are from funnel lines, not the points
fig_HL <- plot(f1) +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom") +
        scale_x_continuous(limits = c(0,350), breaks = c(0,50,100,150,200,250,300,350)) +
        scale_y_continuous(limits = c(0,0.6), breaks = c(0,0.2,0.4,0.6)) +
        labs(x = "", y = "Covariate-adjusted proportion\nreceiving frontline curative RT", subtitle = "A: CHL" ) +
        ggsci::scale_colour_jco(name = "Control limits", labels = c("95%", "99.7%"))

f2 <- funnelModel(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                          imd_2 + imd_3  + imd_4  + imd_5 +
                          chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                          `final_route_Emergency presentation` + final_route_Unknown +
                          `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                          `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                          `distance_km_cat_40+` | Org_of_care ,
                  pointTarget(limits = c(0.05,0.003), standardised = FALSE, normalApprox = TRUE), data = df_analysis_final_DLBCL)

fig_DLBCL <- plot(f2) +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom") +
        scale_x_continuous(limits = c(0,650), breaks = c(0,100,200,300,400,500,600)) +
        scale_y_continuous(limits = c(0,0.6), breaks = c(0,0.2,0.4,0.6,0.8, 1)) +
        labs(x = "Number of patients receiving frontline curative RT", y = "Covariate-adjusted proportion\nreceiving frontline curative RT",
             subtitle = "B: DLBCL" ) +
        ggsci::scale_colour_jco(name = "Control limits", labels = c("95%", "99.7%"))

## Correlation

res_corr_SR <- df_analysis_final_cHL %>%
        group_by(trust.name) %>%
        summarise(observed_HL = sum(RT_binary), expected_HL = sum(expected)) %>%
        mutate(rate_HL = observed_HL/expected_HL) %>%
        full_join(df_analysis_final_DLBCL %>%
                          group_by(trust.name) %>%
                          summarise(observed_DLBCL = sum(RT_binary), expected_DLBCL = sum(expected)) %>%
                          mutate(rate_DLBCL = observed_DLBCL/expected_DLBCL), by = c("trust.name"))

Figure_corr_SR <- ggscatter(res_corr_SR, x = "rate_HL", y = "rate_DLBCL",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Standardised rate: CHL", ylab = "Standardised rate: DLBCL") +
        theme_bw(base_size = 14) +
        labs(subtitle = "C: NHS trust correlation")


ggsave(paste0( results_path, "Figure_2ABC_funnelandcorr.png"),
       patchwork::wrap_plots(patchwork::wrap_plots(fig_HL , fig_DLBCL, ncol = 1) +
                                     patchwork::plot_layout(guides = "collect", axes = "collect_y") & theme(legend.position = 'bottom') ,
                             Figure_corr_SR  , ncol = 1, heights = c(2,1)) ,
       width = 16, height =25, units = "cm")

#' _____________________________________________________________________________
## VPC and MOR for HL and DLBC

#'_____ HL
#'
## Specify model

model_1level_HL <- glm(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                               imd_2 + imd_3  + imd_4  + imd_5 +
                               chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                               `final_route_Emergency presentation` + final_route_Unknown + `age_diag_cat_30-39` + `age_diag_cat_40-49` +
                               `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                               `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                               `distance_km_cat_40+`, family=binomial(link='logit'), data = df_analysis_final_cHL)
summary(model_1level_HL)

model_2level_HL <- glmer(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                                 imd_2 + imd_3  + imd_4  + imd_5 +
                                 chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                                 `final_route_Emergency presentation` + final_route_Unknown + `age_diag_cat_30-39` + `age_diag_cat_40-49` +
                                 `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                                 `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                                 `distance_km_cat_40+` + (1 | trust.code) , family=binomial(link='logit'), data = df_analysis_final_cHL)
summary(model_2level_HL)

res_LiklihoodRatioTest_HL <- lmtest::lrtest(model_1level_HL, model_2level_HL)
res_VPC_HL <- performance::icc(model_2level_HL, by_group = TRUE)
tau2_hospital_HL <- (VarCorr(model_2level_HL)$trust.code[1])^2
MOR_hospital_HL <- exp(sqrt(2*tau2_hospital_HL) * qnorm(0.75))

res_GLM_HL_AME_master <- tidy(model_2level_HL,conf.int=TRUE, exponentiate=TRUE, effects="fixed") %>%
        dplyr::select(term, estimate, conf.low, conf.high, p.value  )

res_GLM_HL_AME <- res_GLM_HL_AME_master %>%
        add_row(term = "MOR", estimate = MOR_hospital_HL, conf.low = NA, conf.high = NA) %>%
        filter(term != "(Intercept)") %>%  #remove row number 1 (The intercept)
        mutate(term = case_when(term == "stage_model_Missing" ~ "Stage: Missing",
                                term == "stage_model_2" ~ "Stage: II",
                                term == "`stage_model_3/4`" ~ "Stage: III/IV",
                                term == "MOR" ~ "MOR (NHS Trust)",
                                term == "imd_2" ~ "IMD: 2nd quintile",
                                term == "imd_3" ~ "IMD: 3rd quintile",
                                term == "imd_4" ~ "IMD: 4th quintile",
                                term == "imd_5" ~ "IMD: 5th quintile",
                                term == "gender_Female" ~ "Sex: Female",
                                term == "distance_km" ~ "Distance to RT centre",
                                term == "`age_diag_cat_30-39`" ~ "Age category: 30-39",
                                term == "`age_diag_cat_40-49`" ~ "Age category: 40-49",
                                term == "`age_diag_cat_50-59`" ~ "Age category: 50-59",
                                term == "`age_diag_cat_60-69`" ~ "Age category: 60-69",
                                term == "`age_diag_cat_70-79`" ~ "Age category: 70-79",
                                term == "`age_diag_cat_80+`" ~ "Age category: 80+",
                                term == "chrl_tot_78_06_1" ~ "CCI score 76-6 months: 1",
                                term == "chrl_tot_78_06_2" ~ "CCI score 76-6 months: 2",
                                term == "chrl_tot_78_06_3" ~ "CCI score 76-6 months: 3",
                                term == "`chrl_tot_78_06_3+`" ~ "CCI score 76-6 months: 4+",
                                term == "`final_route_Emergency presentation`" ~ "Route to diagnosis: Emergency presentation",
                                term == "final_route_Unknown" ~ "Route to diagnosis: Unknown",
                                term == "`distance_km_cat_10-20`" ~ "Distance to nearest RT centre: 10-20km",
                                term == "`distance_km_cat_20-30`" ~ "Distance to nearest RT centre: 20-30km",
                                term == "`distance_km_cat_30-40`" ~ "Distance to nearest RT centre: 30-40km",
                                term == "`distance_km_cat_40+`" ~ "Distance to nearest RT centre: 40+km",
                                .default = term),
               term = factor(term, levels = c("Age category: 30-39",
                                              "Age category: 40-49",
                                              "Age category: 50-59",
                                              "Age category: 60-69",
                                              "Age category: 70-79",
                                              "Age category: 80+",
                                              "Sex: Female","Stage: II",
                                              "Stage: III/IV","Stage: Missing",
                                              "IMD: 2nd quintile",
                                              "IMD: 3rd quintile",
                                              "IMD: 4th quintile",
                                              "IMD: 5th quintile",
                                              "CCI score 76-6 months: 1",
                                              "CCI score 76-6 months: 2",
                                              "CCI score 76-6 months: 3",
                                              "CCI score 76-6 months: 4+",
                                              "Route to diagnosis: Emergency presentation",
                                              "Route to diagnosis: Unknown",
                                             "Distance to nearest RT centre: 10-20km",
                                               "Distance to nearest RT centre: 20-30km",
                                                "Distance to nearest RT centre: 30-40km",
                                             "Distance to nearest RT centre: 40+km",
                                              "MOR (NHS Trust)")))


plot_GLM_HL_AME <- res_GLM_HL_AME %>%
        ggplot(aes(x = forcats::fct_rev(term), y = estimate)) +
        geom_point(shape = 15,
                   size  = 4, width = 0.1,
                   position = "dodge", color="black") +
        geom_errorbar(aes(ymin  = conf.low,
                          ymax  = conf.high),
                      width = 0.2,
                      size  = 0.7,
                      position = "dodge", color="turquoise4") +
        theme(axis.title = element_text(face = "bold")) +
        xlab("") + ylab("Adjusted odds ratios with 95% CI") +
        coord_flip(ylim = c(0, 2.5)) +
        geom_hline(yintercept = 1, color = "red", size = 1) +
        theme_bw() +
        theme(axis.title = element_text(size = 17)) +
        theme(axis.text = element_text(size = 14))


#'_____ dlbcl
#'
model_1level_DLBCL <- glm(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                                  imd_2 + imd_3  + imd_4  + imd_5 +
                                  chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                                  `final_route_Emergency presentation` + final_route_Unknown +
                                  `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                                  `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                                  `distance_km_cat_40+` , family=binomial(link='logit'), data = df_analysis_final_DLBCL)

model_2level_DLBCL <- glmer(RT_binary ~  gender_Female + stage_model_2 + `stage_model_3/4` + stage_model_Missing +
                                    imd_2 + imd_3  + imd_4  + imd_5 +
                                    chrl_tot_78_06_1 +  chrl_tot_78_06_2 + chrl_tot_78_06_3 + `chrl_tot_78_06_3+` +
                                    `final_route_Emergency presentation` + final_route_Unknown +
                                    `age_diag_cat_50-59` + `age_diag_cat_60-69` + `age_diag_cat_70-79` + `age_diag_cat_80+` +
                                    `distance_km_cat_10-20` + `distance_km_cat_20-30`+ `distance_km_cat_30-40` +
                                    `distance_km_cat_40+`  + (1 | trust.code) , family=binomial(link='logit'), data = df_analysis_final_DLBCL)

res_LiklihoodRatioTest_DLBCL <- lmtest::lrtest(model_1level_DLBCL, model_2level_DLBCL)
res_VPC_DLBCL <- performance::icc(model_2level_DLBCL, by_group = TRUE)
tau2_hospital_DLBCL <- (VarCorr(model_2level_DLBCL)$trust.code[1])^2
MOR_hospital_DLBCL <- exp(sqrt(2*tau2_hospital_DLBCL) * qnorm(0.75))

res_GLM_DLBCL_AME_master <- tidy(model_2level_DLBCL,conf.int=TRUE, exponentiate=TRUE, effects="fixed") %>%
        dplyr::select(term, estimate, conf.low, conf.high, p.value )

res_GLM_DLBCL_AME <- res_GLM_DLBCL_AME_master %>%
        add_row(term = "MOR", estimate = MOR_hospital_DLBCL, conf.low = NA, conf.high = NA) %>%
        filter(term != "(Intercept)") %>%  #remove row number 1 (The intercept)
        mutate(term = case_when(term == "stage_model_Missing" ~ "Stage: Missing",
                                term == "stage_model_2" ~ "Stage: II",
                                term == "`stage_model_3/4`" ~ "Stage: III/IV",
                                term == "MOR" ~ "MOR (NHS Trust)",
                                term == "imd_2" ~ "IMD: 2nd quintile",
                                term == "imd_3" ~ "IMD: 3rd quintile",
                                term == "imd_4" ~ "IMD: 4th quintile",
                                term == "imd_5" ~ "IMD: 5th quintile",
                                term == "gender_Female" ~ "Sex: Female",
                                term == "distance_km" ~ "Distance to RT centre",
                                term == "`age_diag_cat_50-59`" ~ "Age category: 50-59",
                                term == "`age_diag_cat_60-69`" ~ "Age category: 60-69",
                                term == "`age_diag_cat_70-79`" ~ "Age category: 70-79",
                                term == "`age_diag_cat_80+`" ~ "Age category: 80+",
                                term == "chrl_tot_78_06_1" ~ "CCI score 76-6 months: 1",
                                term == "chrl_tot_78_06_2" ~ "CCI score 76-6 months: 2",
                                term == "chrl_tot_78_06_3" ~ "CCI score 76-6 months: 3",
                                term == "`chrl_tot_78_06_3+`" ~ "CCI score 76-6 months: 4+",
                                term == "`final_route_Emergency presentation`" ~ "Route to diagnosis: Emergency presentation",
                                term == "final_route_Unknown" ~ "Route to diagnosis: Unknown",
                                term == "`distance_km_cat_10-20`" ~ "Distance to nearest RT centre: 10-20km",
                                term == "`distance_km_cat_20-30`" ~ "Distance to nearest RT centre: 20-30km",
                                term == "`distance_km_cat_30-40`" ~ "Distance to nearest RT centre: 30-40km",
                                term == "`distance_km_cat_40+`" ~ "Distance to nearest RT centre: 40+km",.default = term),
               term = factor(term, levels = c("Age category: 50-59",
                                              "Age category: 60-69",
                                              "Age category: 70-79",
                                              "Age category: 80+",
                                              "Sex: Female","Stage: II",
                                              "Stage: III/IV","Stage: Missing",
                                              "IMD: 2nd quintile",
                                              "IMD: 3rd quintile",
                                              "IMD: 4th quintile",
                                              "IMD: 5th quintile",
                                              "CCI score 76-6 months: 1",
                                              "CCI score 76-6 months: 2",
                                              "CCI score 76-6 months: 3",
                                              "CCI score 76-6 months: 4+",
                                              "Route to diagnosis: Emergency presentation",
                                              "Route to diagnosis: Unknown",
                                              "Distance to nearest RT centre: 10-20km",
                                              "Distance to nearest RT centre: 20-30km",
                                              "Distance to nearest RT centre: 30-40km",
                                              "Distance to nearest RT centre: 40+km",
                                              "MOR (NHS Trust)")))


plot_GLM_DLBCL_AME <- res_GLM_DLBCL_AME %>%
        ggplot(aes(x = forcats::fct_rev(term), y = estimate)) +
        geom_point(shape = 15,
                   size  = 4, width = 0.1,
                   position = "dodge", color="black") +
        geom_errorbar(aes(ymin  = conf.low,
                          ymax  = conf.high),
                      width = 0.2,
                      size  = 0.7,
                      position = "dodge", color="turquoise4") +
        theme(axis.title = element_text(face = "bold")) +
        xlab("") + ylab("Adjusted odds ratios with 95% CI") +
        coord_flip(ylim = c(0, 2.5)) +
        geom_hline(yintercept = 1, color = "red", size = 1) +
        theme_bw() +
        theme(axis.title = element_text(size = 17)) +
        theme(axis.text = element_text(size = 14))

#'_____ dlbcl and HL AME ploted together

plot_AME_both <- rbind( res_GLM_DLBCL_AME %>% mutate(Type = "DLBCL"),
       res_GLM_HL_AME %>% mutate(Type = "cHL")) %>%
        mutate(term = factor(term, levels = c("Age category: 30-39",
                                              "Age category: 40-49",
                                              "Age category: 50-59",
                                              "Age category: 60-69",
                                              "Age category: 70-79",
                                              "Age category: 80+",
                                              "Sex: Female","Stage: II",
                                              "Stage: III/IV","Stage: Missing",
                                              "IMD: 2nd quintile",
                                              "IMD: 3rd quintile",
                                              "IMD: 4th quintile",
                                              "IMD: 5th quintile",
                                              "CCI score 76-6 months: 1",
                                              "CCI score 76-6 months: 2",
                                              "CCI score 76-6 months: 3",
                                              "CCI score 76-6 months: 4+",
                                              "Route to diagnosis: Emergency presentation",
                                              "Route to diagnosis: Unknown",
                                              "Distance to nearest RT centre: 10-20km",
                                              "Distance to nearest RT centre: 20-30km",
                                              "Distance to nearest RT centre: 30-40km",
                                              "Distance to nearest RT centre: 40+km",
                                              "MOR (NHS Trust)"))) %>%
        ggplot(aes(x = forcats::fct_rev(term), y = estimate, group = Type, colour= Type )) +
        geom_point(shape = 15,
                   size  = 4, width = 0.1,
                   position = position_dodge(width = 0.5)) +
        geom_errorbar(aes(ymin  = conf.low,
                          ymax  = conf.high),
                      width = 0.2,
                      size  = 0.7,
                      position = position_dodge(width = 0.5)) +
        theme(axis.title = element_text(face = "bold")) +
        xlab("") + ylab("Adjusted odds ratios with 95% CI") +
        coord_flip(ylim = c(0, 2.5)) +
        geom_hline(yintercept = 1, color = "red", size = 1) +
        theme_bw() +
        theme(axis.title = element_text(size = 17)) +
        theme(axis.text = element_text(size = 12)) +
        ggsci::scale_colour_jco() +
        theme(legend.title = element_blank(),
              legend.position = "bottom")

ggsave(paste0( results_path, "Figure_3_GLM_AME_both.png"),
               plot = plot_AME_both, width = 19, height =19, units = "cm", device='png', dpi=600)

#' _____________________________________________________________________________
## Comparison of terms across the two different models (supplement?)
res_AMEs_all <- data.frame(term  = res_logistic_HL_AME$term) %>%
        # HL LR results
        left_join(res_logistic_HL_AME, by =c("term")  ) %>%
        rename(`Odds ratio` = estimate,
               `p value` = p.value,
               `2.5%` = conf.low,
               `97.5%` = conf.high) %>%
        mutate(across(c(`Odds ratio`,`2.5%`, `97.5%`), \(x) round(x, 2))) %>%
        mutate(across(c(`p value`), \(x) round(x, 3))) %>%
        mutate(log_HL_OR_CI_p = paste0(`Odds ratio`, " (",`2.5%`, " to ", `97.5%`, ", ", `p value`, ")")) %>%
        dplyr::select(term, log_HL_OR_CI_p) %>%
        # HL Mixed effects results
        left_join(res_GLM_HL_AME_master, by =c("term")  ) %>%
        rename(`Odds ratio` = estimate,
               `p value` = p.value,
               `2.5%` = conf.low,
               `97.5%` = conf.high) %>%
        mutate(across(c(`Odds ratio`,`2.5%`, `97.5%`), \(x) round(x, 2))) %>%
        mutate(across(c(`p value`), \(x) round(x, 3))) %>%
        mutate(mix_HL_OR_CI_p = paste0(`Odds ratio`, " (",`2.5%`, " to ", `97.5%`, ", ", `p value`, ")")) %>%
        dplyr::select(term, log_HL_OR_CI_p, mix_HL_OR_CI_p) %>%
        # DLBCL LR results
        left_join(res_logistic_DLBCL_AME, by =c("term")  ) %>%
        rename(`Odds ratio` = estimate,
               `p value` = p.value,
               `2.5%` = conf.low,
               `97.5%` = conf.high) %>%
        mutate(across(c(`Odds ratio`,`2.5%`, `97.5%`), \(x) round(x, 2))) %>%
        mutate(across(c(`p value`), \(x) round(x, 3))) %>%
        mutate(log_DLBCL_OR_CI_p = paste0(`Odds ratio`, " (",`2.5%`, " to ", `97.5%`, ", ", `p value`, ")")) %>%
        dplyr::select(term, log_HL_OR_CI_p, mix_HL_OR_CI_p, log_DLBCL_OR_CI_p) %>%
        # DLBCL Mixed effects results
        left_join(res_GLM_DLBCL_AME_master, by =c("term")  ) %>%
        rename(`Odds ratio` = estimate,
               `p value` = p.value,
               `2.5%` = conf.low,
               `97.5%` = conf.high) %>%
        mutate(across(c(`Odds ratio`,`2.5%`, `97.5%`), \(x) round(x, 2))) %>%
        mutate(across(c(`p value`), \(x) round(x, 3))) %>%
        mutate(mix_DLBCL_OR_CI_p = paste0(`Odds ratio`, " (",`2.5%`, " to ", `97.5%`, ", ", `p value`, ")")) %>%
        dplyr::select(term, log_HL_OR_CI_p, mix_HL_OR_CI_p, log_DLBCL_OR_CI_p, mix_DLBCL_OR_CI_p)
