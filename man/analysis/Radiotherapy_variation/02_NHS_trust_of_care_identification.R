#' _____________________________________________________________________________
#' _____________________________________________________________________________

#### Identifying NHS trust

#' _____________________________________________________________________________
#' _____________________________________________________________________________

library(tidyverse)
library(readstata13)


#' _____________________________________________________________________________
# ## load data

raw_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/Analysis/Data/Raw/2023Update/"

analysed_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/Analysis/Data/Analysis/"

package_data_path <- "K:/CTSUCancer/Projects/PHE_Lymphoma/PHELymphRpackage/man/analysis/"


load(paste0(raw_data_path, "trt_reg.Rdata"))
df_patients <- readstata13::read.dta13(file = paste0(analysed_data_path, "Ptinfo_wide_2024_03_13.dta"))
load(paste0(raw_data_path, "RTDS1.RData"))
load(paste0(raw_data_path, "RTDS2.RData"))
load(paste0(raw_data_path, "RTDS1.RData"))
load(paste0(raw_data_path, "tumour_reg.RData"))

df_tum_short <- tumour_reg %>% dplyr:::select(pseudo_patientid, pseudo_tumourid, diag_trust, diag_trust_name)

df_RTDS <- RTDS1 %>%
        distinct()

df_HL_RTDS <- read.csv(paste0(package_data_path, 'Radiotherapy_variation/',  "HL_RTDS_intent.csv"),
                       na.strings=c("","-","   ","     ")) %>%
        dplyr::select(-X) %>%
        mutate(RT_intent_3 = NA, .before = apptdate_1 ) %>%
        mutate(apptdate_3  = NA, .before = orgcodeprovider_1  ) %>%
        mutate(orgcodeprovider_3 = NA)

df_DLBCL_RTDS <- read.csv(paste0(package_data_path, 'Radiotherapy_variation/',  "DLBCL_RTDS_intent.csv"),
                       na.strings=c("","-","   ","     ")) %>%
        dplyr::select(-X)


df_HL_DLBCL_RTDS <- rbind(df_HL_RTDS, df_DLBCL_RTDS)

#' _____________________________________________________________________________
# ## Create patients with HL and DLBCL dataframe - this is the analysis set

df_RT_analysis <- df_patients %>%
        mutate(pseudo_patientid = as.character(pseudo_patientid)) %>%
        filter(datdiag >= as.Date("2014-01-01")) %>% ## 258647 attendances
        filter(diagcat7 == "CHL" |  diagcat7 == "DLBCL" ) %>%
        left_join(df_HL_DLBCL_RTDS, by = c("pseudo_patientid" = "pseudo_patientid")) %>%
        mutate(pseudo_tumourid = as.character(pseudo_tumourid)) %>%
        left_join(df_tum_short, by = c("pseudo_tumourid"))

df_trt_reg_orgs <- trt_reg %>%
        dplyr::select(pseudo_patientid, pseudo_tumourid, eventdesc, trust_code, trust_name) %>%
        group_by(pseudo_patientid) %>%
        distinct(trust_name, .keep_all=TRUE) %>%
        mutate(event_no=1:n()) %>%
        ungroup() %>%
        filter(event_no < 4) %>%
        pivot_wider(names_from = event_no, values_from = c(eventdesc, trust_code, trust_name)) %>%
        mutate(pseudo_tumourid = as.character(pseudo_tumourid))

df_RT_analysis <- df_RT_analysis %>%
        left_join(df_trt_reg_orgs, by = c("pseudo_tumourid" = "pseudo_tumourid"))

write.csv(df_RT_analysis , file=paste0(package_data_path, 'Radiotherapy_variation/',"df_RT_analysis", '.csv'), row.names = FALSE )


#' _____________________________________________________________________________
#' Load in HES CWT and DIDS

#' Add HES inpatient organisation code for identification

#' As its such a big file run it once, save and comment out
#'
# load(paste0(raw_data_path, "HESinpt.RData"))
# df_HES_inp_orgs_wide <- HESinpt %>%
#         arrange(pseudo_patientid) %>%
#         select(pseudo_patientid,procode3) %>%
#         filter(pseudo_patientid %in% df_RT_analysis$pseudo_patientid.x) %>%
#         group_by(pseudo_patientid) %>%
#         mutate(num.uniq.trusts = n_distinct(procode3)) %>%
#         ungroup() %>%
#         group_by(pseudo_patientid, procode3) %>%
#         mutate(total.num.epps.at.trusts = n()) %>%
#         distinct() %>%
#         ungroup() %>%
#         group_by(pseudo_patientid, num.uniq.trusts) %>%
#         arrange(desc(total.num.epps.at.trusts), .by_group = TRUE) %>%
#         mutate(trust.num = 1:n()) %>%
#         ungroup() %>%
#         group_by(pseudo_patientid) %>%
#         select(-num.uniq.trusts) %>%
#         pivot_wider(names_from = trust.num, values_from = c(procode3, total.num.epps.at.trusts))

#' Save processed version
# write.csv(df_HES_inp_orgs_wide , file=paste0(package_data_path, 'Radiotherapy_variation/', "df_HESInp_orgcodes", '.csv'), row.names = FALSE )

df_RT_analysis <- read.csv(file=paste0(package_data_path, 'Radiotherapy_variation/',"df_RT_analysis.csv"))
df_HES_inp_orgs <- read.csv(file=paste0(package_data_path, 'Radiotherapy_variation/', "df_HESInp_orgcodes.csv"))

#' Join
df_RT_HES_orgs <- df_RT_analysis %>%
        dplyr:::select(pseudo_patientid.x, diag_trust, diag_trust_name, trust_code_1, trust_name_1,
                       orgcodeprovider_1, eventdesc_1:trust_name_3 ) %>%
        left_join(df_HES_inp_orgs, by = c("pseudo_patientid.x" = "pseudo_patientid")) %>%
        mutate(Org_of_care = "Not assigned yet")

#' ________
#' CWT

load(paste0(raw_data_path, "CWT.Rdata"))
df_CWT <- CWT %>%
        group_by(pseudo_patientid) %>%
        slice(1) %>%
        ungroup %>%
        dplyr:::select(pseudo_patientid, org_dec_to_treat) %>%
        mutate(org_dec_to_treat = substr(org_dec_to_treat,1,3))
#' ________
#' DIFD
load(paste0(raw_data_path, "DIDS.Rdata"))
df_DIDs <- DIDS %>%
        dplyr:::select(pseudo_patientid, referring_org_code) %>%
        distinct() %>%
        group_by(pseudo_patientid) %>%
        mutate(event_no=1:n()) %>%
        ungroup() %>%
        filter(event_no < 5) %>%
        pivot_wider(names_from = event_no, names_prefix = "DID" , values_from = c(referring_org_code ))


# Join
df_org_CWT_DIDs <- left_join(df_RT_HES_orgs, df_CWT , by = c("pseudo_patientid.x" = "pseudo_patientid")) %>%
        left_join(df_DIDs, by = c("pseudo_patientid.x" = "pseudo_patientid"))


#' _____________________________________________________________________________
#' Assignment of 'Trust of care' First with CWT (org_dec_to_treat) then cancer reg (diag_trust)
df_org_CWT_DIDs <- df_org_CWT_DIDs %>%
        mutate(across(where(is.character), \(x) na_if(x,""))) %>%
        mutate(Org_of_care = case_when(!is.na(org_dec_to_treat) ~ as.character(org_dec_to_treat), #
                                       !is.na(diag_trust) ~ as.character(diag_trust),
                                     .default =  Org_of_care))


#' Check whether they generally match?
#' org_dec_to_treat = CWT, diag_trust = cancer reg
#' ~72% agreement
df_orgs_check <- df_org_CWT_DIDs %>%
        filter(!is.na(org_dec_to_treat) & !is.na(diag_trust)) %>%
        dplyr::select(pseudo_patientid.x  , org_dec_to_treat, diag_trust) %>%
        mutate(same = org_dec_to_treat==diag_trust) %>%
        summarise(n= n(),
                  n_true = sum(same)) %>%
        mutate(agreement = (n_true/n) *100 )

#
# # Check how many left
# templist <- df_org_CWT_DIDs %>%
#         mutate(Org_of_care = as.factor(Org_of_care )) %>%
#         group_by(Org_of_care) %>%
#         summarise(total = n())

## Filter only those who do not have CWT or DID trust code of diagnosis - only 204 of 39979 patients
df_orgs <- df_org_CWT_DIDs %>%
        filter(Org_of_care =="Not assigned yet")


#' _____________________________________________________________________________
## Best guess based on HES, CAS and DIDs

#' _____________________________________________________________________________
#### Decision fork No.No.3: Does patient have HES codes
No.No.Yes <- df_orgs %>%
        filter(!is.na(total.num.epps.at.trusts_1))

No.No.No <- df_orgs %>%
        filter(is.na(total.num.epps.at.trusts_1))


#' _____________________________________________________________________________
#### Decision fork No.No.Yes.4: Is there a most frequent place? if so, assign

Done.No.No.Yes.Yes <- No.No.Yes %>%
        filter(total.num.epps.at.trusts_1 > total.num.epps.at.trusts_2 | is.na(total.num.epps.at.trusts_2)) %>%  ##
        mutate(Org_of_care = procode3_1)

No.No.Yes.No <- subset(No.No.Yes, !(pseudo_patientid.x %in% Done.No.No.Yes.Yes$pseudo_patientid.x)) ##


#' _____________________________________________________________________________
#### Decision fork No.No.Yes.No.5: CAS code of first event? No = assign HES proceddure code

Done.No.No.Yes.No.Yes <- No.No.Yes.No %>%
        filter(!is.na(trust_code_1)) %>% ##
        mutate(Org_of_care = trust_code_1)

Done.No.No.Yes.No.No <- No.No.Yes.No %>%
        filter(is.na(trust_code_1)) %>% ##
        mutate(Org_of_care = procode3_1)

#'

#' _____________________________________________________________________________
#### Decision fork No.No.No.4: CAS codes?

Done.No.No.No.Yes <- No.No.No %>%
        filter(!is.na(trust_code_1)) %>% ##
        mutate(Org_of_care = trust_code_1)

No.No.No.No <- No.No.No %>%
        filter(is.na(trust_code_1)) ## 24 Exclude from DIDs a lot of private codes)

#' _____________________________________________________________________________
## Commbine dataframes
library(reshape2)

df_list <- mget(ls(pattern = "Done"))
df_complete_orgs <- do.call("rbind", df_list)
rownames(df_complete_orgs) <- NULL


## 24 of 39979 excluded as trust could not be assigned
df_complete_orgs <- rbind.data.frame(df_complete_orgs, df_org_CWT_DIDs %>%
                                             filter(Org_of_care != "Not assigned yet"))

df_complete_orgs <- df_complete_orgs %>%
        dplyr:::select(pseudo_patientid.x, Org_of_care)

write.csv(df_complete_orgs , file=paste0(package_data_path, 'Radiotherapy_variation/',"Patient_NHS_Trust_of_care", '.csv'), row.names = FALSE )

