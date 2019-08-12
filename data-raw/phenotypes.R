library(tidyverse)
options(java.parameters = "-Xmx7g")
library(XLConnect)
library(rlang)

# Working Phenotypes, most phenotypes should be in here
wb_pheno <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                   paste0("GWAS DOE Switchgrass, ",
                                          "Consolidated Data"),
                                   paste0("DOE_GWAS_2019_Consolidated ",
                                          "Pheno Data_Working.xlsx")))
lst_pheno <- readWorksheet(wb_pheno, sheet = getSheets(wb_pheno))
phenos_planting1 <- lst_pheno$`GWAS 2019 Greenup Data`
phenos_planting2 <- lst_pheno$`GWAS 2019 Core Phenotype Data`
phenos_key <- lst_pheno$`Column Key and Notes`

load(file.path("C:", "Users", "ahm543", "OneDrive", "Juenger Projects",
               "SwitchgrassGWAS", "Package", "switchgrassGWAS", "data",
               "metadata.rda"))
load(file.path("C:", "Users", "ahm543", "OneDrive", "Juenger Projects",
               "SwitchgrassGWAS", "Package", "switchgrassGWAS", "data",
               "Taxa.rda"))

# Individual site phenotypes
wb_CLMB2018 <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                      paste0("GWAS DOE Switchgrass, ",
                                             "Consolidated Data"),
                                      "Individual Site Files",
                                      "CLMB_GWAS_2018_Rust Data_Final.xlsx"))
CLMB2018 <- readWorksheet(wb_CLMB2018, sheet = getSheets(wb_CLMB2018))

KBSM2018 <- read_csv(file.path("C:", "Users", "ahm543", "Dropbox",
                               paste0("GWAS DOE Switchgrass, ",
                                      "Consolidated Data"),
                               "Individual Site Files",
                               "KBSM_GWAS_2018_Rust Data_Final.csv"),
                     col_names = TRUE)
KING2018 <- read_csv(file = file.path("..", "..", "..",
                                      "Metadata_and_Phenotypes",
                                      "Armyworm damage KING 2018.csv"),
                     col_names = TRUE)
TMPL2018 <- read_csv(file = file.path("..", "..", "..",
                                      "Metadata_and_Phenotypes",
                                      "Armyworm damage TMPL 2018.csv"),
                     col_names = TRUE)

wb_KING2019 <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                      paste0("GWAS DOE Switchgrass, ",
                                             "Consolidated Data"),
                                      "Individual Site Files",
                                      "KING_GWAS_2019_Chlorosis Data.xlsx"))
lst_KING2019 <- readWorksheet(wb_KING2019, sheet = getSheets(wb_KING2019))
KING2019 <- as_tibble(lst_KING2019$`KING GWAS 2019 Chlorosis`)
key_KING19 <- as_tibble(lst_KING2019$`COLUMN KEY`)

wb_PKLE2018 <- loadWorkbook(file.path("..", "..", "..",
                                      "Metadata_and_Phenotypes",
                                      "PKLE_GWAS_2018_Master_Data_Final.xlsx"))
lst_PKLE2018 = readWorksheet(wb_PKLE2018, sheet = getSheets(wb_PKLE2018))
PKLE_2018 <- as_tibble(lst_PKLE2018$`PKLE_2018_MASTER DATA`)
key_PKLE18 <- as_tibble(lst_PKLE2018$`PKLE_2018_COLUMN KEY and NOTES`)
# phenotype key; this should help make sense of the phenotypic data.
cultural_PKLE18 <- as_tibble(lst_PKLE2018$`PKLE_2018_GWAS_CULTURAL DATA`)
# who did what, when & why; could help make sense of the phenotypic data.

wb_DEAD2018 <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                      paste0("GWAS DOE Switchgrass, ",
                                             "Consolidated Data"),
                                      paste0("GWAS_2019_Mortality, Weak and ",
                                             "Not P. Virgatum Plants.xlsx")))
lst_DEAD2018 <- readWorksheet(wb_DEAD2018, sheet = getSheets(wb_DEAD2018))
DEAD2018 <- as_tibble(lst_DEAD2018$`GWAS 2019 Dead Only`)
LOST2018 <- as_tibble(lst_DEAD2018$`GWAS 2019 Dead, Weak, Not PV`)

LOST2018 <- LOST2018 %>%
  mutate(DEAD_2018 = case_when(
    NOTES %in% c("dead", "Dead", "Dead?", "DEAD") ~ 1,
    NOTES %in% c("1 tiller", "2 tillers", "3 tillers", "4 tillers",
                 "5 tillers", "Very Weak", "small crown", "Partial crown death",
                 "small crown/2 tillers", "small crown/3 tillers",
                 "small crown/4 tillers", "small crown/5 tillers",
                 "small crown/7 tillers", "small crown/8 tillers") ~ 1,
    NOTES %in% c("Not Virgatum", "switchgrass?", "Not Virgatum, Replant") ~ na_dbl,
    TRUE ~ 0
  ))

all_plants_phenos <- phenos_planting1 %>%
  left_join(PKLE_2018, by = c("SITE", "PLOT_GL", "PLOT_LC", "ACC", "INDV", "PLANT_ID")) %>%
  left_join(KBSM2018, by = c("PLOT_GL")) %>%
  left_join(TMPL2018, by = c("SITE", "PLOT_GL", "PLOT_LC")) %>%
  mutate(DAM = DAM/10) %>%
  left_join(KING2018, by = c("SITE", "PLOT_GL", "PLOT_LC")) %>%
  left_join(CLMB2018, by = "PLOT_GL") %>%
  left_join(KING2019, by = "PLOT_GL") %>%
  left_join(LOST2018, by = c("SITE", "PLOT_GL", "PLOT_LC", "PLANT_ID",
                             "PLANT_.ID_GL")) %>%
  mutate(PLANT_ID = ifelse(ACC == "AP13",
                           "AP13",
                           PLANT_ID),
         DEAD_2018 = ifelse(is.na(DEAD_2018), 0, DEAD_2018))

# Need to figure out the best way to join phenos_planting2 here

pvdiv_784g <- metadata %>%
  left_join(all_plants_phenos, by = "PLANT_ID") %>%
  filter(!is.na(SITE)) %>%
  dplyr::select(-(COLLECTION_TYPE:COLLECTION_METHOD))


# Ensure all of the columns have the correct data type - mostly, numeric.

pvdiv_784g <- pvdiv_784g %>%
  mutate(TC_EOS_2018 = as.numeric(TC_EOS_2018),
         SRV_2018 = case_when(
           SRV == "Y" ~ TRUE,
           SRV == "N" ~ FALSE,
           TRUE ~ na_lgl
         ),
         FRZ_DMG = as.numeric(FRZ_DMG),
         LEN = as.numeric(LEN),
         RUST_242 = as.numeric(RUST_242),
         RUST_214 = as.numeric(RUST_214),
         RUST_267 = as.numeric(RUST_267),
         RUST = coalesce(RUST_214, RUST_242),
         RUST = coalesce(RUST, RUST_267),
         Latitude = as.numeric(Latitude),
         Elevation = as.numeric(Elevation),
         LMA = as.numeric(LMA),
         RUST = as.numeric(RUST),
         DAM.x = as.numeric(DAM.x),
         DAM.y = as.numeric(DAM.y),
         DAM = coalesce(DAM.x, DAM.y),
         SLA = as.numeric(SLA),
         MASS_DRY = as.numeric(MASS_DRY),
         AREA = as.numeric(AREA),
         LAM = as.numeric(LAM),
         MID = as.numeric(MID),
         Longitude = as.numeric(Longitude),
         GWAS_CT = as.numeric(GWAS_CT),
         VLAM = LAM*(AREA - (MID*LEN)),
         VMID = pi*LEN*(MID/2)^2,
         ALAM = AREA - (MID*LEN),
         AMID = MID*LEN,
         LVA_GEOM = (VLAM/ALAM) + (VMID/AMID),
         LD_GEOM = LMA/LVA_GEOM,
         LVA_MEAN = (LAM + MID)/2,
         LD_MEAN = LMA/LVA_MEAN,
         LMA_SQD = LMA^2,
         TC_EOS_SQ = sqrt(TC_EOS_2018),
         LAM_SQ = sqrt(LAM),
         MID_SQ = sqrt(MID),
         SLA_RCP = 1/SLA,
         AREA_SQ = sqrt(AREA),
         MASS_SQ = sqrt(MASS_DRY),
         GR1 = as.numeric(GR1),
         GR50 = as.numeric(GR50),
         GR100 = as.numeric(GR100),
         #FL1 = as.numeric(FL1),
         #FL50 = as.numeric(FL50),
         #FL100 = as.numeric(FL100),
         #TC_EOS = as.numeric(TC_EOS),
         #LD_EOS = as.numeric(LD_EOS),
         #HT_PAN_EOS = as.numeric(HT_PAN_EOS),
         #BIOMASS = as.numeric(BIOMASS),
         CHLR_101 = as.numeric(CHLR_101),
         CHLR_121 = as.numeric(CHLR_121),
         CHLR_136 = as.numeric(CHLR_136),
         CHLR_JASON_138 = as.numeric(CHLR_JASON_138),
         CHLR_150. = as.numeric(CHLR_150.),
         SPAD_158 = as.numeric(SPAD_158),
         DEAD_2018 = as.numeric(DEAD_2018)
  )

# Get the mean phenotype value for each PLANT_ID at each SITE as the GWAS value
phenotypes <- pvdiv_784g %>%
  group_by(PLANT_ID, SITE) %>%
  summarise(TC_EOS_2018 = mean(TC_EOS_2018, na.rm = TRUE),
            LEN = mean(LEN, na.rm = TRUE),
            RUST = mean(RUST, na.rm = TRUE),
            Latitude = mean(Latitude, na.rm = TRUE),
            Longitude = mean(Longitude, na.rm = TRUE),
            Elevation = mean(Elevation, na.rm = TRUE),
            LMA = mean(LMA, na.rm = TRUE),
            DAM = mean(DAM, na.rm = TRUE),
            SLA = mean(SLA, na.rm = TRUE),
            MASS_DRY = mean(MASS_DRY, na.rm = TRUE),
            AREA = mean(AREA, na.rm = TRUE),
            LAM = mean(LAM, na.rm = TRUE),
            MID = mean(MID, na.rm = TRUE),
            GWAS_CT = mean(GWAS_CT, na.rm = TRUE),
            VLAM = mean(VLAM, na.rm = TRUE),
            VMID = mean(VMID, na.rm = TRUE),
            ALAM = mean(ALAM, na.rm = TRUE),
            AMID = mean(AMID, na.rm = TRUE),
            LVA_GEOM = mean(LVA_GEOM, na.rm = TRUE),
            LVA_MEAN = mean(LVA_MEAN, na.rm = TRUE),
            LD_MEAN = mean(LD_MEAN, na.rm = TRUE),
            GR1 = mean(GR1, na.rm = TRUE),
            GR50 = mean(GR50, na.rm = TRUE),
            GR100 = mean(GR100, na.rm = TRUE),
            #FL1 = mean(FL1, na.rm = TRUE),
            #FL50 = mean(FL50, na.rm = TRUE),
            #FL100 = mean(FL100, na.rm = TRUE),
            #TC_EOS = mean(TC_EOS, na.rm = TRUE),
            #LD_EOS = mean(LD_EOS, na.rm = TRUE),
            #HT_PAN_EOS = mean(HT_PAN_EOS, na.rm = TRUE),
            #BIOMASS = mean(BIOMASS, na.rm = TRUE),
            CHLR_101 = mean(CHLR_101, na.rm = TRUE),
            CHLR_121 = mean(CHLR_121, na.rm = TRUE),
            CHLR_136 = mean(CHLR_136, na.rm = TRUE),
            CHLR_JASON_138 = mean(CHLR_JASON_138, na.rm = TRUE),
            CHLR_150. = mean(CHLR_150., na.rm = TRUE),
            SPAD_158 = mean(SPAD_158, na.rm = TRUE),
            DEAD_2018 = min(DEAD_2018, na.rm = TRUE)
  ) %>%
  mutate_all( ~ case_when(!is.nan(.x) ~ .x)) %>% # if the column is not NA, keep the value; else replace with NA as nothing is provided as another choice in case_when.
  filter(!is.na(SITE) & !is.na(PLANT_ID))

#NB: No outlier removal was done here.
phenotypes <- phenotypes %>%
  gather(TC_EOS_2018:RUST,LMA:MID,VLAM:DEAD_2018, key = "Phenotype", value = value) %>%
  filter(!is.na(value)) %>%
  # saveRDS(file = "Cleaned pvdiv Phenotypes from 2018 for 784g 2019-05-16.rds")
  unite(col = "SITE_PHE", SITE, Phenotype) %>%
  spread(key = SITE_PHE, value = value) %>%
  right_join(Taxa)

#  saveRDS(phenotypes, file = "upland-vs-lowland/phenotypes.rds")
#  GR10 <- phenotypes %>%
#    dplyr::select(PLANT_ID, ends_with("GR1"), ends_with("GR50"), ends_with("GR100"))
#  save(GR10, file = "greenup_2019_30_phenotypes.rda")
#  load("C:/Users/alice/OneDrive/Juenger Projects/SwitchgrassGWAS/Package/switchgrassGWAS/data/phenotypes.rda")
phenotypes <- phenotypes %>%
  dplyr::select(PLANT_ID, GWAS_CT, BRKG_TC_EOS_2018, BRKG_DEAD_2018)

usethis::use_data(phenotypes, compress = "gzip", overwrite = TRUE)
