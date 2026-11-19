#' _____________________________________________________________________________
#' _____________________________________________________________________________

#### HL

#' _____________________________________________________________________________
#' _____________________________________________________________________________

library(tidyverse)
library(readstata13)

#' _____________________________________________________________________________
#' ______ load data

raw_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/Analysis/Data/Raw/2023Update/"

analysed_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/Analysis/Data/Analysis/"

package_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/PHELymphRpackage/man/analysis/"

# For Ruhi: The trt_reg, and RTDs datasets are the original NCRAS datasets. The Ptinfo_wide_2024_03_13.dta file is a dataset created by another member of David Cutter's team
# which is an almalgomation of the NCRAS tumour, patient, and treatment datasets
load(paste0(raw_data_path, "trt_reg.Rdata"))
df_patients <- readstata13::read.dta13(file = paste0(analysed_data_path, "Ptinfo_wide_2024_03_13.dta"))
load(paste0(raw_data_path, "RTDS1.RData"))
load(paste0(raw_data_path, "RTDS2.RData"))

# remove duplicate attendences
df_RTDS <- RTDS1 %>%
        distinct()

# Join patient info - only actually need diagnosis date, diagcat7 which a 7 category grouping of lymphoma types
df_RTDS_Master <- df_RTDS %>%
        left_join(df_patients %>% mutate(pseudo_patientid = as.character(pseudo_patientid)), by = c("pseudo_patientid" = "pseudo_patientid")) %>%
        mutate(across(contains("dat"), as.Date, format = "%d/%m/%Y", origin="1970-01-01")) %>%
        mutate(datediff_diag_trt = difftime(treatmentstartdate, datdiag, units = "days")) %>%
        filter(grepl("P", pseudo_patientid))

#' _____________________________________________________________________________
#' ______ HL datafrane
df_HL <- df_RTDS_Master %>% ##
        filter(treatmentstartdate >= as.Date("2014-01-01")) %>% ##
        filter(diagcat7 == "CHL" ) %>%
        mutate(RT_intent = "Not assigned yet") %>%
        group_by(pseudo_patientid, pseudo_rtepisodeid) %>%
        mutate(num.app = n()) %>%
        mutate(num.uniq.prescdose = n_distinct(rtprescribeddose)) %>%
        mutate(sum.uniq.prescdose = sum(unique(rtprescribeddose), na.rm = TRUE)) %>%
        mutate(num.uniq.prescfractions = n_distinct(prescribedfractions)) %>%
        mutate(sum.uniq.prescfractions = sum(unique(prescribedfractions), na.rm = TRUE)) %>%
        mutate(date_diff_trt_death = difftime(datdth, treatmentstartdate, units = "days")) %>%
        ungroup()

#' _____________________________________________________________________________
#### Decision fork 1: Treatment more than a year after diagnosis?
### only looking at variation in first line RT,
## everyone with RT a year afer diagnosis simply coded as "Later lines or palliative"

Done.Yes <- df_HL %>%
        filter(datediff_diag_trt > 365) %>%
        mutate(RT_intent = "Later lines or palliative")

No <- df_HL %>%
        filter(datediff_diag_trt <= 365)


#' _____________________________________________________________________________
#### Diagnosis less  than a year ###############################################
#' _____________________________________________________________________________
#### Decision fork No.2: Curative regimen?

No.Yes <- df_HL %>%
        filter(datediff_diag_trt <= 365) %>%
        filter((between(round(rtprescribeddose, 0), 25, 40) & between(round(prescribedfractions, 0), 5, 20)) |
                       (between(round(rtprescribeddose, 0), 18, 24) & between(round(prescribedfractions, 0), 6, 15)) ) %>%
        mutate(RT_intent = "1L Curative")


No.No <- df_HL %>%
        filter(datediff_diag_trt <= 365) %>%
        filter(!((between(round(rtprescribeddose, 0), 25, 40) & between(round(prescribedfractions, 0), 5, 20)) |
                         (between(round(rtprescribeddose, 0), 18, 24) & between(round(prescribedfractions, 0), 6, 15))) |
                       is.na(rtprescribeddose) | is.na(prescribedfractions))

#' _____________________________________________________________________________
#### Decision fork No.Yes.3: use sum dose and fractions

Done.No.Yes.Yes <- No.Yes %>%
        filter(grepl(c("C8"),radiotherapydiagnosisicd)) %>%
        mutate(RT_intent = "1L Curative")

No.Yes.No <- No.Yes %>%
        filter(!(grepl(c("C8"),radiotherapydiagnosisicd)))

#' _____________________________________________________________________________
#### Decision fork No.Yes.No.4: use sum dose and fractions

Done.No.Yes.No.yes <- No.Yes.No %>%
        filter((substr(radiotherapydiagnosisicd, start = 1, stop = 2) == substr(earlier_diag_01, start = 1, stop = 2) |
                       substr(radiotherapydiagnosisicd, start = 1, stop = 2) == substr(later_diag_01, start = 1, stop = 2)) &
                       rtprescribeddose != 30)

Done.No.Yes.No.No <- No.Yes.No %>%
        filter(!((substr(radiotherapydiagnosisicd, start = 1, stop = 2) == substr(earlier_diag_01, start = 1, stop = 2) |
                          substr(radiotherapydiagnosisicd, start = 1, stop = 2) == substr(later_diag_01, start = 1, stop = 2)) &
                         rtprescribeddose != 30))


#' _____________________________________________________________________________
#### Decision fork No.No.3: use sum dose and fractions

Done.No.No.Yes <- No.No %>%
        filter((between(round(sum.uniq.prescdose, 0), 25, 40) & between(round(sum.uniq.prescfractions, 0), 5, 20)) |
                       (between(round(sum.uniq.prescdose, 0), 18, 24) & between(round(sum.uniq.prescfractions, 0), 6, 15))) %>%
        mutate(RT_intent = "1L Curative")

No.No.No <- No.No %>%
        filter(!((between(round(sum.uniq.prescdose, 0), 25, 40) & between(round(sum.uniq.prescfractions, 0), 5, 20)) |
                         (between(round(sum.uniq.prescdose, 0 ), 18, 24) & between(round(sum.uniq.prescfractions, 0), 6, 15))) |
                       is.na(sum.uniq.prescdose) | is.na(sum.uniq.prescfractions))


#' _____________________________________________________________________________
#### Decision fork No.No.No.4: Ealier or later cancer?

Done.No.No.No.Yes <- No.No.No %>%
        filter((!is.na(earlier_diagdate_01)| !is.na(later_diagdate_01)) & !grepl("C81",radiotherapydiagnosisicd ))  %>%
        mutate(RT_intent = "Not related to index lymphoma")

No.No.No.No <- No.No.No %>%
        filter(!((!is.na(earlier_diagdate_01)| !is.na(later_diagdate_01)) & !grepl("C81",radiotherapydiagnosisicd )))

#' _____________________________________________________________________________
#### Decision fork No.No.No.NO.5:  Palliative regimen?

Done.No.No.No.No.Yes <- No.No.No.No  %>%
        filter((between(round(rtprescribeddose, 0 ), 20, 20) & between(round(prescribedfractions, 0), 5, 5)) |
                       (between(round(sum.uniq.prescdose, 0), 20, 20) & between(round(sum.uniq.prescfractions, 0), 5, 5)) |
                       (between(round(rtprescribeddose, 0), 0, 17) ) |
                       (between(round(sum.uniq.prescdose, 0), 0, 17))) %>%
        mutate(RT_intent = "Later lines or palliative")



No.No.No.No.No <- No.No.No.No %>%
        filter(!((between(round(rtprescribeddose, 0 ), 20, 20) & between(round(prescribedfractions, 0), 5, 5)) |
                         (between(round(sum.uniq.prescdose, 0), 20, 20) & between(round(sum.uniq.prescfractions, 0), 5, 5)) |
                         (between(round(rtprescribeddose, 0), 0, 17) ) |
                         (between(round(sum.uniq.prescdose, 0), 0, 17))) |
                       is.na(rtprescribeddose) | is.na(prescribedfractions))


#' _____________________________________________________________________________
#### Decision fork No.No.No.NO.No.6:  For a different cancer

Done.No.No.No.No.No.Yes <- No.No.No.No.No  %>%
        filter((!grepl(c("C81"),radiotherapydiagnosisicd) | (!is.na(earlier_diagdate_01)| !is.na(later_diagdate_01))) &
                       rtprescribeddose != 30) %>%
        mutate(RT_intent = "Not related to index lymphoma")


Done.No.No.No.No.No.No <- No.No.No.No.No %>%
        filter(!((!grepl(c("C81"),radiotherapydiagnosisicd) | (!is.na(earlier_diagdate_01)| !is.na(later_diagdate_01))) &
                       rtprescribeddose != 30)) %>%
        mutate(RT_intent = "1L Curative")



#' _____________________________________________________________________________
#### Final HL datset ###############################################
#' _____________________________________________________________________________
## Check those with multiple prescriptions/doses/fractions per episode and then
## make df_HL_RTDs one line per patient

library(reshape2)

df_list <- mget(ls(pattern = "Done"))
df_HL_intent <- do.call("rbind", df_list) ## Adds to the right number of attendances


df_HL_intent <- df_HL_intent %>%
        mutate_if(is.character,as.factor)

# First check each attendance for a patient episode has same RT intent. Make two dataframes. First a dataframe
# for those who had different intents. Remove the rows with the wrong intent. Then a second for
# those with the same intent. Make it one attendances per patient episode than bind dataframes back together

df_check_HL_1intent <- df_HL_intent %>%
        relocate(RT_intent, .before = 1) %>%
        group_by(pseudo_patientid, pseudo_rtepisodeid) %>%
        mutate(sum.uniq.intent = n_distinct(RT_intent)) %>%
        filter(sum.uniq.intent > 1) %>%
        ungroup() %>%
        group_by(pseudo_patientid, pseudo_rtepisodeid, RT_intent) %>%
        filter(row_number() == 1) %>%
        #group_by(pseudo_patientid) %>% summarise(n = n()) # check
        arrange(pseudo_patientid) %>%
        # Rules based on looking at data - not logic rules hence weird C7 condition
        filter((RT_intent == "Not related to index lymphoma" & !grepl(c("C8"), radiotherapydiagnosisicd)) |
                       (RT_intent == "1L Curative" & grepl(c("C8|C7"), radiotherapydiagnosisicd))) %>%
        #group_by(pseudo_patientid) %>% summarise(n = n()) # check
        select(pseudo_patientid, pseudo_rtepisodeid, RT_intent,apptdate, orgcodeprovider)


df_check_HL_2intent <- df_HL_intent %>%
        group_by(pseudo_patientid, pseudo_rtepisodeid) %>%
        mutate(sum.uniq.intent = n_distinct(RT_intent)) %>%
        filter(sum.uniq.intent == 1) %>%
        filter(row_number() == 1)%>%
        select(pseudo_patientid, pseudo_rtepisodeid, RT_intent,apptdate, orgcodeprovider)

df_HL_intent_small <- rbind(df_check_HL_1intent, df_check_HL_2intent)

# Now see how many have multiple episodes. And for those with multiple episodes whether they
# have the same intent.

# Those with one episode
df_check_HL_1intent <- df_HL_intent_small %>%
        group_by(pseudo_patientid) %>%
        mutate(sum.epp = n_distinct(pseudo_rtepisodeid),
               sum.intent = n_distinct(RT_intent)) %>%
        filter(sum.epp == 1)

# For those with multiple episodes but with same intent, take only the first episode
df_check_HL_2intent <- df_HL_intent_small %>%
        group_by(pseudo_patientid) %>%
        mutate(sum.epp = n_distinct(pseudo_rtepisodeid),
               sum.intent = n_distinct(RT_intent)) %>%
        filter(sum.epp > 1) %>%
        filter(sum.intent == 1) %>%
        filter(row_number() == 1)

# For those with multiple intent, remove nonsensical ones -
# only ePID = 651554.12460.11306 - palliative before curative. Then one row per intent.
# E.g there may be three episodes, curative curative palliative. Make into curative palliative


df_check_HL_3intent <- df_HL_intent_small %>%
        group_by(pseudo_patientid) %>%
        mutate(sum.epp = n_distinct(pseudo_rtepisodeid),
               sum.intent = n_distinct(RT_intent)) %>%
        filter(sum.epp > 1) %>%
        filter(sum.intent > 1) %>%
        arrange(pseudo_patientid, apptdate) %>%
        group_by(pseudo_patientid) %>%
        mutate(RT_intent2 = lag(RT_intent)) %>%
        filter(!(RT_intent == "1L Curative" & RT_intent2 == "Later lines or palliative") | is.na(RT_intent2)) %>%
        group_by(pseudo_patientid, RT_intent) %>%
        filter(row_number() == 1)

# Join back together and pivot so one row per patient. Check it match correct number
#
nrow(df_HL %>% group_by(pseudo_patientid) %>%
        filter(row_number() == 1)) ## 3333

df_HL_intent_small_2 <- rbind(df_check_HL_1intent, df_check_HL_2intent, df_check_HL_3intent)
df_HL_intent_small_2 <- df_HL_intent_small_2 %>%
        select(pseudo_patientid, pseudo_rtepisodeid, RT_intent, apptdate, orgcodeprovider) %>%
        group_by(pseudo_patientid) %>%
        mutate(epp.num=1:n()) %>%
        select(-pseudo_rtepisodeid) %>%
        pivot_wider(names_from = epp.num, values_from = c("RT_intent", "apptdate", "orgcodeprovider"))

# 3333 matches!

write.csv(df_HL_intent_small_2, file=paste0(package_data_path, 'Radiotherapy_variation/',  'HL_RTDS_intent', '.csv'))


