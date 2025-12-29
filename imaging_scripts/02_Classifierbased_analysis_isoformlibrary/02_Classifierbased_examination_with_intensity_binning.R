### R script to perform classifier-based analysis of protein isoform pairs with intensity-based binning ###
# Note, in the script below, three different binning approaches are possible. We tried all three and they resulted in similar results.
# We ultimately used APPROACH 1 and 3 in the study as they allowed for the most balanced control of intensity levels between isoforms within a given pair within the bins. 
# Note, In the final study, Approach 3 was relabelled Approach 2, and approach 2 (as specified below) is not included

#packages
library(tidyverse)
library(caret)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(broom)
library(readxl)
library(lsr) #Contains cohensD() function
library(nnet) #Feedforward neural network
library(UpSetR) #Package for doing upset plots
library(ggrepel) #Pacakge for making repel text

### Functions for Loading filtered transfection datasets AND creating collected summary ###
#Function for loading in TRANSFECTED CELLS and SUMMARY datasets for each SET.
Load_transfected_cells_and_summary <- function(SET){
  setwd(paste0('INSERT DIRECORY'))
 
  Transfected_cell_files <- list.files() %>% 
    str_subset(pattern = '_TRANSFECTED_CELLS')
  Transfection_summary_files <- list.files() %>% 
    str_subset(pattern = '_TRANSFECTION_SUMMARY')
  
  load(paste0(str_subset(Transfected_cell_files, pattern = paste0(SET,'_REP_1_TRANSFECTED_CELLS'))))
  TRANSFECTED_CELLS_REP_1 <- Transfected_cellbodies
  rm(Transfected_cellbodies)
  
  load(paste0(str_subset(Transfected_cell_files, pattern = paste0(SET,'_REP_2_TRANSFECTED_CELLS'))))
  TRANSFECTED_CELLS_REP_2 <- Transfected_cellbodies
  rm(Transfected_cellbodies)
  
  load(paste0(str_subset(Transfected_cell_files, pattern = paste0(SET,'_REP_3_TRANSFECTED_CELLS'))))
  TRANSFECTED_CELLS_REP_3 <- Transfected_cellbodies
  rm(Transfected_cellbodies)
  
  Transfected_cells <- bind_rows(TRANSFECTED_CELLS_REP_1,
                                 TRANSFECTED_CELLS_REP_2,
                                 TRANSFECTED_CELLS_REP_3)
  
  load(paste0(str_subset(Transfection_summary_files, pattern = paste0(SET,'_REP_1_TRANSFECTION_SUMMARY'))))
  TRANSFECTION_SUMMARY_REP_1 <- Transfection_summary
  rm(Transfection_summary)
  
  load(paste0(str_subset(Transfection_summary_files, pattern = paste0(SET,'_REP_2_TRANSFECTION_SUMMARY'))))
  TRANSFECTION_SUMMARY_REP_2 <- Transfection_summary
  rm(Transfection_summary)
  
  load(paste0(str_subset(Transfection_summary_files, pattern = paste0(SET,'_REP_3_TRANSFECTION_SUMMARY'))))
  TRANSFECTION_SUMMARY_REP_3 <- Transfection_summary
  rm(Transfection_summary)
  
  Transfection_summary <- bind_rows(TRANSFECTION_SUMMARY_REP_1,
                                    TRANSFECTION_SUMMARY_REP_2,
                                    TRANSFECTION_SUMMARY_REP_3)
  return(list(Transfected_cells, Transfection_summary))
}
#Function for collected summary of REPEAT summaries
COLLECTED_SUMMARY <- function(DATASET){
  coll_summ <- DATASET %>%  
    ungroup() %>% 
    group_by(Experiment_Set, Metadata_Well) %>% 
    summarise(Transfected = sum(Transfected_TRUE),
              n_cells = sum(n_cells),
              Transfection_rate = round(Transfected/n_cells * 100, digits = 1)) %>%
    dplyr::mutate('Transfected_REP_1' = DATASET$Transfected_TRUE[DATASET$Experiment_Rep == 'REP_1'],
                  'Transfected_REP_2' = DATASET$Transfected_TRUE[DATASET$Experiment_Rep == 'REP_2'],
                  'Transfected_REP_3' = DATASET$Transfected_TRUE[DATASET$Experiment_Rep == 'REP_3'],
                  'n_cells_REP_1' = DATASET$n_cells[DATASET$Experiment_Rep == 'REP_1'],
                  'n_cells_REP_2' = DATASET$n_cells[DATASET$Experiment_Rep == 'REP_2'],
                  'n_cells_REP_3' = DATASET$n_cells[DATASET$Experiment_Rep == 'REP_3'],
                  'Transfection_rate_REP_1' = round(Transfected_REP_1/n_cells_REP_1 * 100, digits = 1),
                  'Transfection_rate_REP_2' = round(Transfected_REP_2/n_cells_REP_2 * 100, digits = 1),
                  'Transfection_rate_REP_3' = round(Transfected_REP_3/n_cells_REP_3 * 100, digits = 1))
  return(coll_summ)
}

### Special case for SET 1 / A as they mostly contain the same isoforms
#AGGREGATING SET 1 AND SET A
SET_1_FOR_FUSION <- Load_transfected_cells_and_summary(SET = 'SET_1')[[1]] %>% 
  dplyr::filter(!c(SET_PAIR %in% c('SET_1_PAIR_06')))
#
SET_A_FOR_FUSION <- Load_transfected_cells_and_summary(SET = 'SET_A')[[1]]
TRANSFECTED_SET_1Aagg <- bind_rows(SET_1_FOR_FUSION,
                                   SET_A_FOR_FUSION)


#
REP_SUMMARY_SET_1_for_fusion <- Load_transfected_cells_and_summary(SET = 'SET_1')[[2]] %>% 
  dplyr::filter(!c(SET_PAIR %in% c('SET_1_PAIR_06')))
#
REP_SUMMARY_SET_A_for_fusion <- Load_transfected_cells_and_summary(SET = 'SET_A')[[2]]

COLLECTED_SUMMARY_SET_1_fusion <- COLLECTED_SUMMARY(REP_SUMMARY_SET_1_for_fusion)
COLLECTED_SUMMARY_SET_A_fusion <- COLLECTED_SUMMARY(REP_SUMMARY_SET_A_for_fusion)
COLLECTED_SUMMARY_SET_1Aagg <- bind_rows(COLLECTED_SUMMARY_SET_1_fusion,
                                         COLLECTED_SUMMARY_SET_A_fusion)
#In this case, PAIR06 WILL BE VAV2 (AS PAIR06 IS REMOVE FROM SET_1. YOU THEN NEED TO RENAME THE SET_PAIR VARIABLE)

### Loading transfected sets for rest.

#SET A
TRANSFECTED_SET_A <- Load_transfected_cells_and_summary(SET = 'SET_A')[[1]]
REP_SUMMARY_SET_A <- Load_transfected_cells_and_summary(SET = 'SET_A')[[2]]
#SET B
TRANSFECTED_SET_B <- Load_transfected_cells_and_summary(SET = 'SET_B')[[1]]
REP_SUMMARY_SET_B <- Load_transfected_cells_and_summary(SET = 'SET_B')[[2]]
#SET C
TRANSFECTED_SET_C <- Load_transfected_cells_and_summary(SET = 'SET_C')[[1]]
REP_SUMMARY_SET_C <- Load_transfected_cells_and_summary(SET = 'SET_C')[[2]]
#SET D
TRANSFECTED_SET_D <- Load_transfected_cells_and_summary(SET = 'SET_D')[[1]] %>% 
  dplyr::filter(!c(SET_PAIR == 'SET_D_PAIR_NegativeGolgin97'))
REP_SUMMARY_SET_D <- Load_transfected_cells_and_summary(SET = 'SET_D')[[2]] %>% 
  dplyr::filter(!c(SET_PAIR == 'SET_D_PAIR_NegativeGolgin97'))
#SET E
TRANSFECTED_SET_E <- Load_transfected_cells_and_summary(SET = 'SET_E')[[1]] %>% 
  dplyr::mutate(SET_PAIR = if_else(SET_PAIR == 'SET_E_PAIR_NegativeGolgin97','Set_D_PAIR_NEGATIVE', SET_PAIR))

REP_SUMMARY_SET_E <- Load_transfected_cells_and_summary(SET = 'SET_E')[[2]] %>% 
  dplyr::mutate(SET_PAIR = if_else(SET_PAIR == 'SET_E_PAIR_NegativeGolgin97','Set_D_PAIR_NEGATIVE', SET_PAIR))


COLLECTED_SUMMARY_SET_A <- COLLECTED_SUMMARY(REP_SUMMARY_SET_A)
COLLECTED_SUMMARY_SET_B <- COLLECTED_SUMMARY(REP_SUMMARY_SET_B)
COLLECTED_SUMMARY_SET_C <- COLLECTED_SUMMARY(REP_SUMMARY_SET_C)
COLLECTED_SUMMARY_SET_D <- COLLECTED_SUMMARY(REP_SUMMARY_SET_D)
COLLECTED_SUMMARY_SET_E <- COLLECTED_SUMMARY(REP_SUMMARY_SET_E)
#
COLLECTED_SUMMARY_SET_1 <- COLLECTED_SUMMARY(REP_SUMMARY_SET_1_for_fusion)

#####
### Removing variables to reduce storage requirements. Currently removes all non-mCherry variables. Note, can be made to keep variables if REDUCE_FEATURES = FALSE
#Function to remove variables. ONLY_TRANSFECTED takes FALSE/TRUE as input. If TRUE, filters dataset for transfected cells
remove_non_mcherry_vars <- function(TRANSFECTED_SET, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = FALSE){
  if(REDUCE_FEATURES == TRUE){
    Reduced_dataset <- TRANSFECTED_SET %>% 
      dplyr::select(ImageNumber, ObjectNumber, OBJECT_ID, Metadata_Well, CONTROL,Experiment_Set, Experiment_Rep, SET_PAIR, ISOFORM,TRANSFECTED, 
                    starts_with('Intensity_') & ends_with('_MCHERRY_BACKCORR'),
                    starts_with('Correlation_'),
                    starts_with('Granularity_'), starts_with('RadialDistribution'), starts_with('Texture_'), starts_with('AreaShape_')) %>% 
      dplyr::select(!c(starts_with('Correlation'))) %>% 
      na.omit()
  }
  if(REDUCE_FEATURES == FALSE){
    Reduced_dataset <- TRANSFECTED_SET %>% 
      dplyr::select(ImageNumber, ObjectNumber, OBJECT_ID, Metadata_Well, CONTROL,Experiment_Set, Experiment_Rep, SET_PAIR, ISOFORM,TRANSFECTED, everything()) %>% 
      dplyr::select(!c(Experiment_SetRep, ends_with('_BACKGROUND'), Number_Object_Number, Parent_Segmentation)) %>% 
      na.omit()
  }
  if(ONLY_TRANSFECTED == TRUE){
    Reduced_dataset <- Reduced_dataset %>% 
      dplyr::filter(TRANSFECTED == TRUE)
  }
  if(ONLY_NEGTRANSFECTED == TRUE){
    Reduced_dataset <- Reduced_dataset %>% 
      dplyr::filter(TRANSFECTED == FALSE)
  }
  return(Reduced_dataset)
}
#
Reduced_transfected_1 <- remove_non_mcherry_vars(SET_1_FOR_FUSION, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE)
Reduced_transfected_A <- remove_non_mcherry_vars(TRANSFECTED_SET_A, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE)
Reduced_transfected_B <- remove_non_mcherry_vars(TRANSFECTED_SET_B, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE)
Reduced_transfected_C <- remove_non_mcherry_vars(TRANSFECTED_SET_C, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE) %>% 
  dplyr::filter(!c(SET_PAIR %in% c('SET_C_PAIR_01','SET_C_PAIR_03'))) 
Reduced_transfected_D <- remove_non_mcherry_vars(TRANSFECTED_SET_D, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE)
Reduced_transfected_1Aagg <- remove_non_mcherry_vars(TRANSFECTED_SET_1Aagg, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE)
Reduced_transfected_E <- remove_non_mcherry_vars(TRANSFECTED_SET_E, ONLY_TRANSFECTED = TRUE, REDUCE_FEATURES = TRUE)

#
Reduced_negtrans_1 <- remove_non_mcherry_vars(SET_1_FOR_FUSION, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE)
Reduced_negtrans_A <- remove_non_mcherry_vars(TRANSFECTED_SET_A, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE)
Reduced_negtrans_B <- remove_non_mcherry_vars(TRANSFECTED_SET_B, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE)
Reduced_negtrans_C <- remove_non_mcherry_vars(TRANSFECTED_SET_C, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE) %>% 
  dplyr::filter(!c(SET_PAIR %in% c('SET_C_PAIR_01','SET_C_PAIR_03')))
Reduced_negtrans_D <- remove_non_mcherry_vars(TRANSFECTED_SET_D, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE)
Reduced_negtrans_1Aagg <- remove_non_mcherry_vars(TRANSFECTED_SET_1Aagg, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE)
Reduced_negtrans_E <- remove_non_mcherry_vars(TRANSFECTED_SET_E, ONLY_TRANSFECTED = FALSE, REDUCE_FEATURES = TRUE, ONLY_NEGTRANSFECTED = TRUE)

#Alternative pairings of isoforms
#For Set A and Set D/E
 #Defined here as may be used in intensity calculations
#####
#Alternative comparisons (can be used for both 1, A, and 1A)
alternative_comparisons_set_1A <- function(DATASET, SET_NAME){
  alt_pairClasp_altEx <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'DEx' | SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'DEx') %>% 
    dplyr::mutate(
      ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'DEx' ~ 'DEx',
                          SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'DEx' ~ 'PanEx'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_03DExVSPanEx'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pairClasp_bothex <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'DEx' | SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'WT') %>% 
    dplyr::mutate(ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'DEx' ~ 'DEx',
                                      SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'WT' ~ 'Both'),
                  SET_PAIR = paste0(SET_NAME,'_PAIR_03bothEx'),
                  SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pairClasp_exvsex <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'WT' | SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'DEx') %>% 
    dplyr::mutate(
      ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'WT' ~ 'Ex1',
                          SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'DEx' ~ 'PanEx'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_03ExVsEx'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pairClasp_ex1vsboth <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'WT' | SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'WT') %>% 
    dplyr::mutate(
      ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_03') & ISOFORM == 'WT' ~ 'Ex1',
                          SET_PAIR == paste0(SET_NAME,'_PAIR_04') & ISOFORM == 'WT' ~ 'Both'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_03Ex1VsBoth'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  Transfected_only_without_outliers_1A_ALTERNATIVE_COMPARISONS <- DATASET %>% 
    dplyr::filter(SET_PAIR %in% c('SET_1_PAIR_MCHERRY','SET_1_PAIR_RAB5A', 'SET_A_PAIR_MCHERRY','SET_A_PAIR_RAB5A','SET_1A_PAIR_MCHERRY','SET_1A_PAIR_RAB5A')) %>% 
    bind_rows(., alt_pairClasp_altEx, alt_pairClasp_bothex, alt_pairClasp_exvsex,alt_pairClasp_ex1vsboth) %>% 
    na.omit()
  
  return(Transfected_only_without_outliers_1A_ALTERNATIVE_COMPARISONS)
}

#For SET D (ADDITIONAL COMPARISONS)
#Single comparisons (for App and Vav2) in Set D / E
alternative_comparisons_set_D <- function(DATASET, SET_NAME){ # DATASET is the transfected cells, SET_NAME is either SET_D or SET_E (Other sets don't have specific pairs)
  alt_pair01_wt57 <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_01') & ISOFORM %in% c('WT', 'DEx57nt')) %>% 
    dplyr::mutate(
      ISOFORM = case_when(ISOFORM == 'WT' ~ 'WT',
                          ISOFORM == 'DEx57nt' ~ 'DEx'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_01wtVS57'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pair01_wt168 <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_01') & ISOFORM %in% c('WT', 'DEx168nt')) %>% 
    dplyr::mutate(
      ISOFORM = case_when(ISOFORM == 'WT' ~ 'WT',
                          ISOFORM == 'DEx168nt' ~ 'DEx'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_01wtVS168'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pair01_dexdex <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_01') & ISOFORM %in% c('DEx57nt', 'DEx168nt')) %>% 
    dplyr::mutate(
      ISOFORM = case_when(ISOFORM == 'DEx57nt' ~ 'DEx57nt',
                          ISOFORM == 'DEx168nt' ~ 'DEx168nt'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_01dex57VS168'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  #
  alt_pair08altex <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'DEx' | SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'DEx') %>% 
    dplyr::mutate(
      ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'DEx' ~ 'DEx',
                          SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'DEx' ~ 'AltEx'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_08DExVSaltEx'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  
  alt_pair08bothex <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'DEx' | SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'WT') %>% 
    dplyr::mutate(ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'DEx' ~ 'DEx',
                                      SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'WT' ~ 'Both'),
                  SET_PAIR = paste0(SET_NAME,'_PAIR_08bothEx'),
                  SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pair08exvsex <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'WT' | SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'DEx') %>% 
    dplyr::mutate(
      ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'WT' ~ 'Ex1',
                          SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'DEx' ~ 'Ex2'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_08ExVsEx'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  alt_pair08ex1vsboth <- DATASET %>% 
    dplyr::filter(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'WT' | SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'WT') %>% 
    dplyr::mutate(
      ISOFORM = case_when(SET_PAIR == paste0(SET_NAME,'_PAIR_08') & ISOFORM == 'WT' ~ 'Ex1',
                          SET_PAIR == paste0(SET_NAME,'_PAIR_09') & ISOFORM == 'WT' ~ 'Both'),
      SET_PAIR = paste0(SET_NAME,'_PAIR_08Ex1VsBoth'),
      SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  #Creating dataset with alterantive comparisons (and controls)
  Transfected_only_without_outliers_D_ALTERNATIVE_COMPARISONS <- DATASET %>% 
    dplyr::filter(SET_PAIR %in% c('SET_D_PAIR_MCHERRY','SET_D_PAIR_RAB5A', 'SET_E_PAIR_MCHERRY','SET_E_PAIR_RAB5A')) %>% 
    bind_rows(., alt_pair01_wt57, alt_pair01_wt168, alt_pair01_dexdex, alt_pair08altex, alt_pair08bothex, alt_pair08exvsex, alt_pair08ex1vsboth) %>% 
    na.omit()
  
  return(Transfected_only_without_outliers_D_ALTERNATIVE_COMPARISONS)
}

### Creating a summary of the relative median intensity differences between isoforms (for plotting to diagnose classifier)
#####

# PRE-FUSE 1 and A to make counts / intesity stuff easier for color plotting / diagnostics (not used for actual classifier input although probably would work fine)
#####
# Postrans
Reduced_transfected_1Aagg_fused <- Reduced_transfected_1Aagg %>% 
  dplyr::mutate(Experiment_Rep = case_when(Experiment_Set == 'SET_1' & Experiment_Rep == 'REP_1' ~ 'REP_1',
                                           Experiment_Set == 'SET_1' & Experiment_Rep == 'REP_2' ~ 'REP_2',
                                           Experiment_Set == 'SET_1' & Experiment_Rep == 'REP_3' ~ 'REP_3',
                                           Experiment_Set == 'SET_A' & Experiment_Rep == 'REP_1' ~ 'REP_4',
                                           Experiment_Set == 'SET_A' & Experiment_Rep == 'REP_2' ~ 'REP_5',
                                           Experiment_Set == 'SET_A' & Experiment_Rep == 'REP_3' ~ 'REP_6')) %>% 
  dplyr::mutate(Experiment_Set = 'SET_1A') %>%
  dplyr::mutate(SET_PAIR = paste0(Experiment_Set,'_',word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))),
                SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM)) %>% 
  dplyr::filter(!c(ISOFORM == 'NEGATIVE')) %>% 
  na.omit()
# Negtrans
Reduced_negtrans_1Aagg_fused <- Reduced_negtrans_1Aagg %>% 
  dplyr::mutate(Experiment_Rep = case_when(Experiment_Set == 'SET_1' & Experiment_Rep == 'REP_1' ~ 'REP_1',
                                           Experiment_Set == 'SET_1' & Experiment_Rep == 'REP_2' ~ 'REP_2',
                                           Experiment_Set == 'SET_1' & Experiment_Rep == 'REP_3' ~ 'REP_3',
                                           Experiment_Set == 'SET_A' & Experiment_Rep == 'REP_1' ~ 'REP_4',
                                           Experiment_Set == 'SET_A' & Experiment_Rep == 'REP_2' ~ 'REP_5',
                                           Experiment_Set == 'SET_A' & Experiment_Rep == 'REP_3' ~ 'REP_6')) %>% 
  dplyr::mutate(Experiment_Set = 'SET_1A') %>%
  dplyr::mutate(SET_PAIR = paste0(Experiment_Set,'_',word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))),
                SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM)) %>% 
  dplyr::filter(!c(ISOFORM == 'NEGATIVE')) %>% 
  na.omit()

#
### Set 1A, alternative comparisons
# PositiveTransfectants
Reduced_transfected_1Aagg_fused_ALTERNATIVE_COMPARISONS <- alternative_comparisons_set_1A(Reduced_transfected_1Aagg_fused, SET_NAME = 'SET_1A')
# NegativeTransfectants
Reduced_negtrans_1Aagg_fused_ALTERNATIVE_COMPARISONS <- alternative_comparisons_set_1A(Reduced_negtrans_1Aagg_fused, SET_NAME = 'SET_1A')


#D alternative comparisons, prefused before intensity filtering stuff
Reduced_transfected_D_alt_comps <- alternative_comparisons_set_D(DATASET = Reduced_transfected_D, SET_NAME = 'SET_D') %>%  #NEEDS THE ALT COMP SET EFUNCTIION ()
  dplyr::mutate(ISOFORM = case_when(SET_PAIR == 'SET_D_PAIR_01wtVS57' & ISOFORM == 'WT' ~ 'WT',
                                    SET_PAIR == 'SET_D_PAIR_01wtVS57' & ISOFORM == 'DEx' ~ 'DEx',
                                    SET_PAIR == 'SET_D_PAIR_01wtVS168' & ISOFORM == 'WT' ~ 'WT',
                                    SET_PAIR == 'SET_D_PAIR_01wtVS168' & ISOFORM == 'DEx' ~ 'DEx',
                                    SET_PAIR == 'SET_D_PAIR_01dex57VS168' & ISOFORM == 'DEx57nt' ~ 'WT',
                                    SET_PAIR == 'SET_D_PAIR_01dex57VS168' & ISOFORM == 'DEx168nt' ~ 'DEx',
                                    ISOFORM == 'WT' ~ 'WT',
                                    ISOFORM == 'DEx' ~ 'DEx',
                                    ISOFORM == 'AltEx' ~ 'WT',
                                    ISOFORM == 'Both' ~ 'WT',
                                    ISOFORM == 'Ex1' ~ 'DEx',
                                    ISOFORM == 'Ex2' ~ 'WT'))
Reduced_negtrans_D_alt_comps <- alternative_comparisons_set_D(DATASET = Reduced_negtrans_D, SET_NAME = 'SET_D') %>%  #NEEDS THE ALT COMP SET EFUNCTIION ()
  dplyr::mutate(ISOFORM = case_when(SET_PAIR == 'SET_D_PAIR_01wtVS57' & ISOFORM == 'WT' ~ 'WT',
                                    SET_PAIR == 'SET_D_PAIR_01wtVS57' & ISOFORM == 'DEx' ~ 'DEx',
                                    SET_PAIR == 'SET_D_PAIR_01wtVS168' & ISOFORM == 'WT' ~ 'WT',
                                    SET_PAIR == 'SET_D_PAIR_01wtVS168' & ISOFORM == 'DEx' ~ 'DEx',
                                    SET_PAIR == 'SET_D_PAIR_01dex57VS168' & ISOFORM == 'DEx57nt' ~ 'WT',
                                    SET_PAIR == 'SET_D_PAIR_01dex57VS168' & ISOFORM == 'DEx168nt' ~ 'DEx',
                                    ISOFORM == 'WT' ~ 'WT',
                                    ISOFORM == 'DEx' ~ 'DEx',
                                    ISOFORM == 'AltEx' ~ 'WT',
                                    ISOFORM == 'Both' ~ 'WT',
                                    ISOFORM == 'Ex1' ~ 'DEx',
                                    ISOFORM == 'Ex2' ~ 'WT'))
#####



### Creating positive control ###
# Only one per set, so no need to do in alt comparisons
# Function
PositiveControl_comparisons <- function(DATASET, SET_NAME){
  poscontrol_comparison <- DATASET %>% 
    dplyr::filter(Metadata_Well == 'A02' | Metadata_Well == 'C01' | Metadata_Well == 'C02') %>% 
    dplyr::mutate(
      #Note, only C02 of thE RAB5A wells are used for the PosCOntrol so only two wells contribute the cells
      ISOFORM = case_when(Metadata_Well == 'A02' ~ 'WT',
                          Metadata_Well == 'C02' ~ 'DEx',
                          Metadata_Well == 'C01' ~ 'WT'),
      SET_PAIR = case_when(Metadata_Well == 'A02' ~ paste0(SET_NAME,'_PAIR_PosControl'),
                           Metadata_Well == 'C02' ~ paste0(SET_NAME,'_PAIR_PosControl'),
                           Metadata_Well == 'C01' ~ paste0(SET_NAME,'_PAIR_ALONE'))) %>% 
    dplyr::filter(Metadata_Well != 'C01')
  return(poscontrol_comparison)
}

# SET 1
Reduced_transfected_1_w_controls <- bind_rows(Reduced_transfected_1,
                                                       PositiveControl_comparisons(DATASET = Reduced_transfected_1, SET_NAME = 'SET_1'))

Reduced_negtrans_1_w_controls <- bind_rows(Reduced_negtrans_1,
                                                          PositiveControl_comparisons(DATASET = Reduced_negtrans_1, SET_NAME = 'SET_1'))  

# SET A
Reduced_transfected_A_w_controls <- bind_rows(Reduced_transfected_A,
                                                       PositiveControl_comparisons(DATASET = Reduced_transfected_A, SET_NAME = 'SET_A'))       
Reduced_negtrans_A_w_controls <- bind_rows(Reduced_negtrans_A,
                                                          PositiveControl_comparisons(DATASET = Reduced_negtrans_A, SET_NAME = 'SET_A'))  

# SET 1A
Reduced_transfected_1A_w_controls <- bind_rows(Reduced_transfected_1Aagg_fused,
                                                        PositiveControl_comparisons(DATASET = Reduced_transfected_1Aagg_fused, SET_NAME = 'SET_1A'))
Reduced_negtrans_1A_w_controls <- bind_rows(Reduced_negtrans_1Aagg_fused,
                                                           PositiveControl_comparisons(DATASET = Reduced_negtrans_1Aagg_fused, SET_NAME = 'SET_1A'))

# SET B
Reduced_transfected_B_w_controls <- bind_rows(Reduced_transfected_B,
                                                       PositiveControl_comparisons(DATASET = Reduced_transfected_B, SET_NAME = 'SET_B'))       
Reduced_negtrans_B_w_controls <- bind_rows(Reduced_negtrans_B,
                                                          PositiveControl_comparisons(DATASET = Reduced_negtrans_B, SET_NAME = 'SET_B'))  
# SET C
Reduced_transfected_C_w_controls <- bind_rows(Reduced_transfected_C,
                                                       PositiveControl_comparisons(DATASET = Reduced_transfected_C, SET_NAME = 'SET_C'))       
Reduced_negtrans_C_w_controls <- bind_rows(Reduced_negtrans_C,
                                                          PositiveControl_comparisons(DATASET = Reduced_negtrans_C, SET_NAME = 'SET_C'))  
# SET D
Reduced_transfected_D_w_controls <- bind_rows(Reduced_transfected_D,
                                                       PositiveControl_comparisons(DATASET = Reduced_transfected_D, SET_NAME = 'SET_D'))       
Reduced_negtrans_D_w_controls <- bind_rows(Reduced_negtrans_D,
                                                          PositiveControl_comparisons(DATASET = Reduced_negtrans_D, SET_NAME = 'SET_D'))  



#####


#Outliers (DEFINED AS BEING AS HAVING (*CURRENTLY DEFINED AS PER ISOFORM PER REP (SO WELL-BASED):
#Intensity_MeanIntensity_MCHERRY_BACKCORR > 1.5 x IQR + 75% quantile Intensity_MCHERRY_BACKCORR)
###
#Note, function also removes NAN columns and zero colummns
Identify_and_remove_outliers <- function(DATASET,
                                         USE_LOWER_MAX_UPPER_HINGE,
                                         ONLY_NEGTRANS = FALSE,
                                         REMOVE_TOP_x_PERCENTILE_NEGATIVE = FALSE,
                                         TOP_x_PERCENTILE = 0.1,
                                         REMOVE_SAME_COLS_AS_TRANSFECTED = FALSE,
                                         TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = NULL,
                                         relative_intensity_threshold_for_binning = 0.2, # Variable to define what % of max intensity defines lower intensity bin
                                         relative_intensity_percentile_bins = 0.5 # Variable to define what percentiles the bins are split into when cells are bined based on percentile intensity rather than x% of max
                                         ){
  TRANSFECTED_ONLY <- DATASET %>% 
    dplyr::filter(TRANSFECTED == TRUE)
  #In the case that only nontransfected cells are in the original dataset, we don't want to remove them by filtering (as we don't have any cells left then)
  if(ONLY_NEGTRANS == TRUE){
    TRANSFECTED_ONLY <- DATASET 
  }
  if(ONLY_NEGTRANS == TRUE & REMOVE_TOP_x_PERCENTILE_NEGATIVE == TRUE){
    #Read in dataset
    TRANSFECTED_ONLY <- DATASET %>% 
      # Group by the variables of interest
      group_by(SET_PAIR, Experiment_Rep, ISOFORM) %>%
      # Arrange within groups by intensity (descending)
      arrange(desc(Intensity_MeanIntensity_MCHERRY_BACKCORR), .by_group = TRUE) %>%
      # Calculate the percentile within each group
      mutate(percentile = row_number()/n()) %>%
      # Keep only the bottom 80%
      filter(percentile > TOP_x_PERCENTILE) %>%
      # Remove the helper column
      select(-percentile) %>%
      # Ungroup to remove grouping structure
      ungroup()
  }
  
  SUMMARY_STATISTICS <- TRANSFECTED_ONLY %>%
    dplyr::mutate('LOG_INTENSITY' = log10(Intensity_MeanIntensity_MCHERRY_BACKCORR)) %>% 
    ungroup() %>% 
    group_by(SET_PAIR, Experiment_Rep, ISOFORM) %>% 
    summarise(n_cells = length(SET_PAIR),
              MEAN_INTENSITY = mean(Intensity_MeanIntensity_MCHERRY_BACKCORR),
              MEDIAN_INTENSITY = mean(Intensity_MedianIntensity_MCHERRY_BACKCORR),
              SD_INTENSITY = sd(Intensity_MeanIntensity_MCHERRY_BACKCORR),
              LOG_MEAN_INTENSITY = mean(LOG_INTENSITY),
              LOG_SD_INTENSITY = sd(LOG_INTENSITY),
              UPPER_HINGE_IQRx1.5 = quantile(Intensity_MeanIntensity_MCHERRY_BACKCORR, probs = 0.75) + 1.5 * (quantile(Intensity_MeanIntensity_MCHERRY_BACKCORR, probs = 0.75) - quantile(Intensity_MeanIntensity_MCHERRY_BACKCORR, probs = 0.25))) %>% 
    dplyr::mutate('SET_PAIR_REP_ISOFORM' = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  #Can be used to remove the skew of very bright cells if present in one isoform but not another (simply allows using the highest UPPER hing value from one isoform as the filter for the other isoform)
  MAX_INTENSITY_CUTOFF_SET_BASED <- SUMMARY_STATISTICS %>%
    ungroup() %>% 
    group_by(SET_PAIR, ISOFORM) %>%
    summarise('MAX_UPPER_HINGE' = max(UPPER_HINGE_IQRx1.5)) %>% 
    ungroup() %>% 
    group_by(SET_PAIR) %>% 
    summarise('LOWER_MAX_UPPER_HINGE' = min(MAX_UPPER_HINGE))
  #
  
  SUMMARY_STATISTICS_OUTLIER_CRITERIUM <- SUMMARY_STATISTICS %>%
    ungroup() %>% 
    dplyr::select(SET_PAIR_REP_ISOFORM,
                  UPPER_HINGE_IQRx1.5)
  
  Transfected_cells_outliers_removed <- TRANSFECTED_ONLY %>%
    dplyr::mutate('SET_PAIR_REP_ISOFORM' = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM),
                  'LOG_INTENSITY' = log10(Intensity_MeanIntensity_MCHERRY_BACKCORR)) %>% 
    dplyr::mutate('Intensity_Edge_mean_norm' = Intensity_MeanIntensityEdge_MCHERRY_BACKCORR / Intensity_MeanIntensity_MCHERRY_BACKCORR,
                  'Intensity_Edge_max_norm' = Intensity_MeanIntensityEdge_MCHERRY_BACKCORR / Intensity_MeanIntensity_MCHERRY_BACKCORR,
                  'Intensity_Max_mean_norm' = Intensity_MaxIntensity_MCHERRY_BACKCORR / Intensity_MeanIntensity_MCHERRY_BACKCORR,
                  'Intensity_MassDisp_mean_norm' = Intensity_MassDisplacement_MCHERRY_BACKCORR / Intensity_MeanIntensity_MCHERRY_BACKCORR,
                  'Intensity_Std_mean_norm' = Intensity_StdIntensity_MCHERRY_BACKCORR / Intensity_MeanIntensity_MCHERRY_BACKCORR,
                  'Intensity_MAD_mean_norm' = Intensity_MADIntensity_MCHERRY_BACKCORR / Intensity_MeanIntensity_MCHERRY_BACKCORR) %>% 
    dplyr::left_join(x = ., y = SUMMARY_STATISTICS_OUTLIER_CRITERIUM, by = c('SET_PAIR_REP_ISOFORM' = 'SET_PAIR_REP_ISOFORM')) %>%
    dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR < UPPER_HINGE_IQRx1.5) %>%
    dplyr::select(c(1:10,'SET_PAIR_REP_ISOFORM','LOG_INTENSITY','UPPER_HINGE_IQRx1.5', everything()))
  
  if(USE_LOWER_MAX_UPPER_HINGE == TRUE){
    Transfected_cells_outliers_removed <- Transfected_cells_outliers_removed %>% 
      dplyr::left_join(x = ., y = MAX_INTENSITY_CUTOFF_SET_BASED, by = c('SET_PAIR' = 'SET_PAIR')) %>% 
      dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR < LOWER_MAX_UPPER_HINGE)
  }

  ### Creating intensity bin cell subsets ###
  if(ONLY_NEGTRANS == FALSE){
    ### Approach 1 ###
    intensity_bin_thresholds_ap1 <- Transfected_cells_outliers_removed %>%
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      dplyr::summarize(lowest_transfected = min(Intensity_MeanIntensity_MCHERRY_BACKCORR),
              max_transfected = max(Intensity_MeanIntensity_MCHERRY_BACKCORR)) %>% 
      dplyr::mutate(lower_bin_upper_threshold = (relative_intensity_threshold_for_binning * (max_transfected - lowest_transfected)) + lowest_transfected) %>% 
      ungroup()
    #
    Transfected_cells_outliers_removed_ap1_lowerbin <- Transfected_cells_outliers_removed %>%
      left_join(x = ., y = intensity_bin_thresholds_ap1, by = join_by(SET_PAIR)) %>% 
      dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR >= lowest_transfected & Intensity_MeanIntensity_MCHERRY_BACKCORR <= lower_bin_upper_threshold) %>% 
      #
      dplyr::select(-all_of(c('lowest_transfected', 'max_transfected', 'lower_bin_upper_threshold')))
    #
    Transfected_cells_outliers_removed_ap1_upperbin <- Transfected_cells_outliers_removed %>%
      left_join(x = ., y = intensity_bin_thresholds_ap1, by = join_by(SET_PAIR)) %>% 
      dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR > lower_bin_upper_threshold & Intensity_MeanIntensity_MCHERRY_BACKCORR <= max_transfected) %>% 
      #
      dplyr::select(-all_of(c('lowest_transfected', 'max_transfected', 'lower_bin_upper_threshold')))
    
    ### Approach 2 ###
    Transfected_cells_outliers_removed_ap2_lowerbin <- Transfected_cells_outliers_removed %>%
      ungroup() %>% 
      group_by(SET_PAIR, ISOFORM) %>%
      # Arrange within groups by intensity (descending)
      arrange(desc(Intensity_MeanIntensity_MCHERRY_BACKCORR), .by_group = TRUE) %>%
      # Calculate the percentile within each group
      mutate(percentile = row_number()/n()) %>%
      # Keep only the bottom x% percentile
      filter(percentile >= relative_intensity_percentile_bins) %>%
      # Remove the helper column
      select(-percentile) %>%
      # Ungroup to remove grouping structure
      ungroup()
    #
    Transfected_cells_outliers_removed_ap2_upperbin <- Transfected_cells_outliers_removed %>%
      ungroup() %>% 
      group_by(SET_PAIR, ISOFORM) %>%
      # Arrange within groups by intensity (descending)
      arrange(desc(Intensity_MeanIntensity_MCHERRY_BACKCORR), .by_group = TRUE) %>%
      # Calculate the percentile within each group
      mutate(percentile = row_number()/n()) %>%
      # Keep only the bottom x% percentile
      filter(percentile < relative_intensity_percentile_bins) %>%
      # Remove the helper column
      select(-percentile) %>%
      # Ungroup to remove grouping structure
      ungroup()
    
    ### Approach 3 ###
    intensity_bin_thresholds_ap3 <- Transfected_cells_outliers_removed %>%
      ungroup() %>% 
      group_by(SET_PAIR, ISOFORM) %>%
      # Arrange within groups by intensity (descending)
      arrange(desc(Intensity_MeanIntensity_MCHERRY_BACKCORR), .by_group = TRUE) %>%
      # Calculate the percentile within each group
      mutate(percentile = row_number()/n()) %>%
      # Keep only the bottom x% percentile
      filter(percentile >= relative_intensity_percentile_bins) %>%
      # 
      slice_head(n = 1) %>% 
      #
      dplyr::select(SET_PAIR, ISOFORM, Intensity_MeanIntensity_MCHERRY_BACKCORR) %>% 
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      summarize(intensity_threshold_ap3 = max(Intensity_MeanIntensity_MCHERRY_BACKCORR)) %>% 
      ungroup()
  
    #
    Transfected_cells_outliers_removed_ap3_lowerbin <- Transfected_cells_outliers_removed %>%
      ungroup() %>%
      left_join(x = ., y = intensity_bin_thresholds_ap3, by = join_by(SET_PAIR)) %>% 
      dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= intensity_threshold_ap3)
    #
    Transfected_cells_outliers_removed_ap3_upperbin <- Transfected_cells_outliers_removed %>%
      ungroup() %>%
      left_join(x = ., y = intensity_bin_thresholds_ap3, by = join_by(SET_PAIR)) %>% 
      dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR > intensity_threshold_ap3)
   }

  #Seems to remove columns with NANs and no variance
  reduce_cols <- function(DATASET){
    FEATURE_START <- 14
    SET_PAIRS <- DATASET %>% 
      dplyr::filter(!c(word(SET_PAIR, 4,4, sep = fixed('_')) == 'NEGATIVE'),
                    !c(word(SET_PAIR, 4,4, sep = fixed('_')) == 'MCHERRY')) %>% 
      dplyr::pull(SET_PAIR) %>% 
      unique()
    CONSTANT_COLS <-c()
    NAN_COLS <- c()
    FEATURES <- colnames(DATASET)
    FEATURE_END <- length(colnames(DATASET))
    
    #
    for(i in SET_PAIRS){
      PAIR_x_df <- DATASET %>% 
        dplyr::filter(SET_PAIR == i)
      
      for(j in FEATURE_START:FEATURE_END){
        constant_value <- length(unique(PAIR_x_df[[j]]))
        NAN_value <- any(is.na(PAIR_x_df[[j]]))
        
        if(constant_value == 1){
          CONSTANT_COLS <- c(CONSTANT_COLS, j)
        }
        if(NAN_value == TRUE){
          NAN_COLS <- c(NAN_COLS, j)
        }
      }
    }
    CONSTANT_FEATURES <- FEATURES[CONSTANT_COLS]
    NAN_FEATURES <- FEATURES[NAN_COLS]
    FEATURES_TO_REMOVE <- c(unique(CONSTANT_FEATURES), unique(NAN_FEATURES)) #REMOVES 8 VARIABLES FOR SET A (MOSTLY Manders Correlation)
    TEXT_CONSTANT <- paste0('The variables with CONSTANT values for at least one PAIR are: ', unique(CONSTANT_FEATURES))
    TEXT_NULL <- paste0('The variables with MISSING values for at least one PAIR are: ', unique(NAN_FEATURES))
    print(TEXT_CONSTANT)
    print(TEXT_NULL)
    
    Transfected_cells_outliers_removed_REDUCED_COLS_from_function <- DATASET %>% 
      dplyr::select(!c(all_of(FEATURES_TO_REMOVE)))
    
    return(Transfected_cells_outliers_removed_REDUCED_COLS_from_function)
  }
  #
  Transfected_cells_outliers_removed_REDUCED_COLS <- reduce_cols(DATASET = Transfected_cells_outliers_removed)

  if(ONLY_NEGTRANS == FALSE){
    Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp <- reduce_cols(DATASET = Transfected_cells_outliers_removed_ap1_lowerbin)
    Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp <- reduce_cols(DATASET = Transfected_cells_outliers_removed_ap1_upperbin)
    Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp <- reduce_cols(DATASET = Transfected_cells_outliers_removed_ap2_lowerbin)
    Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp <- reduce_cols(DATASET = Transfected_cells_outliers_removed_ap2_upperbin)
    Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp <- reduce_cols(DATASET = Transfected_cells_outliers_removed_ap3_lowerbin)
    Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp <- reduce_cols(DATASET = Transfected_cells_outliers_removed_ap3_upperbin)
 
    #
    keep_only_shared_cols_all <- function(DATASET_1, DATASET_2, DATASET_3, DATASET_4, DATASET_5, DATASET_6){
      common_cols_1 <- intersect(colnames(DATASET_1), colnames(DATASET_2))
      common_cols_2 <- intersect(colnames(DATASET_3), colnames(DATASET_4))
      common_cols_3 <- intersect(colnames(DATASET_5), colnames(DATASET_6))
      #
      common_cols_all <- intersect(common_cols_3, intersect(common_cols_1, common_cols_2))
      #
      DATASET_1_with_common_cols <- DATASET_1 %>% 
        dplyr::select(all_of(common_cols_all))
      #
      DATASET_2_with_common_cols <- DATASET_2 %>% 
        dplyr::select(all_of(common_cols_all))
      #
      DATASET_3_with_common_cols <- DATASET_3 %>% 
        dplyr::select(all_of(common_cols_all))
      #
      DATASET_4_with_common_cols <- DATASET_4 %>% 
        dplyr::select(all_of(common_cols_all))
      #
      DATASET_5_with_common_cols <- DATASET_5 %>% 
        dplyr::select(all_of(common_cols_all))
      #
      DATASET_6_with_common_cols <- DATASET_6 %>% 
        dplyr::select(all_of(common_cols_all))
      
      #
      return(list(DATASET_1_with_common_cols, DATASET_2_with_common_cols,
                  DATASET_3_with_common_cols, DATASET_4_with_common_cols,
                  DATASET_5_with_common_cols, DATASET_6_with_common_cols))
    }
    #
    Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS <- keep_only_shared_cols_all(DATASET_1 = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_2 = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_3 = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_4 = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_5 = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_6 = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp)[[1]]
    Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS <- keep_only_shared_cols_all(DATASET_1 = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_2 = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_3 = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_4 = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_5 = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_6 = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp)[[2]]
    #
    Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS <- keep_only_shared_cols_all(DATASET_1 = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_2 = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_3 = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_4 = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_5 = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_6 = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp)[[3]]
    
    Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS <- keep_only_shared_cols_all(DATASET_1 = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_2 = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_3 = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_4 = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_5 = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_6 = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp)[[4]]
    #
    Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS <- keep_only_shared_cols_all(DATASET_1 = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_2 = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_3 = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_4 = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_5 = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_6 = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp)[[5]]
    Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS <- keep_only_shared_cols_all(DATASET_1 = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_2 = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_3 = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_4 = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS_temp,
                                                                                          DATASET_5 = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS_temp,
                                                                                          DATASET_6 = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS_temp)[[6]]
  }
  
  # For generating negative dataset
  if(ONLY_NEGTRANS == TRUE && REMOVE_SAME_COLS_AS_TRANSFECTED == TRUE){
    Transfected_cells_outliers_removed_REDUCED_COLS <- Transfected_cells_outliers_removed %>%
     dplyr::select(c(colnames(TRANSFECTED_DATASET_FOR_COLS_TO_KEEP)))
  }
  
  # May have to update if we want to calculate imbalance for the subsetted datasets, but don't think it is necessary atm
  SUMMARY_STATISTICS_OUTLIERS_REMOVED <- Transfected_cells_outliers_removed %>%
    ungroup() %>% 
    group_by(SET_PAIR, Experiment_Rep, ISOFORM) %>% 
    summarise(n_cells = length(SET_PAIR),
              MEAN_INTENSITY = mean(Intensity_MeanIntensity_MCHERRY_BACKCORR),
              MEDIAN_INTENSITY = mean(Intensity_MedianIntensity_MCHERRY_BACKCORR),
              SD_INTENSITY = sd(Intensity_MeanIntensity_MCHERRY_BACKCORR),
              LOG_MEAN_INTENSITY = mean(LOG_INTENSITY),
              LOG_SD_INTENSITY = sd(LOG_INTENSITY),
              UPPER_HINGE_IQRx1.5 = quantile(Intensity_MeanIntensity_MCHERRY_BACKCORR, probs = 0.75) + 1.5 * (quantile(Intensity_MeanIntensity_MCHERRY_BACKCORR, probs = 0.75) - quantile(Intensity_MeanIntensity_MCHERRY_BACKCORR, probs = 0.25))) %>% 
    dplyr::mutate('SET_PAIR_REP_ISOFORM' = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
  
  #
  if(ONLY_NEGTRANS == TRUE){
    output_list <- list('full_dataset' = Transfected_cells_outliers_removed_REDUCED_COLS, 'summary_stats' = SUMMARY_STATISTICS, 'summary_stats_without_outliers' = SUMMARY_STATISTICS_OUTLIERS_REMOVED)
  }
  if(ONLY_NEGTRANS == FALSE){
    output_list <- list('full_dataset' = Transfected_cells_outliers_removed_REDUCED_COLS, 'summary_stats' = SUMMARY_STATISTICS, 'summary_stats_without_outliers' = SUMMARY_STATISTICS_OUTLIERS_REMOVED,
                        'ap1_lowerbin' = Transfected_cells_outliers_removed_ap1_lowerbin_REDUCED_COLS, 'ap1_upperbin' = Transfected_cells_outliers_removed_ap1_upperbin_REDUCED_COLS,
                        'ap2_lowerbin' = Transfected_cells_outliers_removed_ap2_lowerbin_REDUCED_COLS, 'ap2_upperbin' = Transfected_cells_outliers_removed_ap2_upperbin_REDUCED_COLS,
                        'ap3_lowerbin' = Transfected_cells_outliers_removed_ap3_lowerbin_REDUCED_COLS, 'ap3_upperbin' = Transfected_cells_outliers_removed_ap3_upperbin_REDUCED_COLS)
  }
  return(output_list)
}

#


# Transfected cells for each isoform pair and controls with binning
## Set 1
Transfected_only_without_outliers_1 <- Identify_and_remove_outliers(Reduced_transfected_1_w_controls, 
                                                                    USE_LOWER_MAX_UPPER_HINGE = TRUE)#
# Note, granularity 14 is removed as SET 1 is the only example where transfected retains more variables than nontrans and this way gives us the same cols
#Approach 1
Transfected_for_class_1_ap1_lowerbin <- Transfected_only_without_outliers_1$ap1_lowerbin %>% 
  dplyr::select(-Granularity_14_MCHERRY_BACKCORR)
Transfected_for_class_1_ap1_upperbin <- Transfected_only_without_outliers_1$ap1_upperbin %>% 
  dplyr::select(-Granularity_14_MCHERRY_BACKCORR)
#Approach 2
Transfected_for_class_1_ap2_lowerbin <- Transfected_only_without_outliers_1$ap2_lowerbin %>% 
  dplyr::select(-Granularity_14_MCHERRY_BACKCORR)
Transfected_for_class_1_ap2_upperbin <- Transfected_only_without_outliers_1$ap2_upperbin %>% 
  dplyr::select(-Granularity_14_MCHERRY_BACKCORR)
#Approach 3
Transfected_for_class_1_ap3_lowerbin <- Transfected_only_without_outliers_1$ap3_lowerbin %>% 
  dplyr::select(-Granularity_14_MCHERRY_BACKCORR)
Transfected_for_class_1_ap3_upperbin <- Transfected_only_without_outliers_1$ap3_upperbin %>% 
  dplyr::select(-Granularity_14_MCHERRY_BACKCORR)


## Set A
Transfected_only_without_outliers_A <- Identify_and_remove_outliers(Reduced_transfected_A_w_controls, 
                                                                    USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_A_ap1_lowerbin <- Transfected_only_without_outliers_A$ap1_lowerbin
Transfected_for_class_A_ap1_upperbin <- Transfected_only_without_outliers_A$ap1_upperbin
#Approach 2
Transfected_for_class_A_ap2_lowerbin <- Transfected_only_without_outliers_A$ap2_lowerbin
Transfected_for_class_A_ap2_upperbin <- Transfected_only_without_outliers_A$ap2_upperbin
#Approach 3
Transfected_for_class_A_ap3_lowerbin <- Transfected_only_without_outliers_A$ap3_lowerbin
Transfected_for_class_A_ap3_upperbin <- Transfected_only_without_outliers_A$ap3_upperbin


## Set 1A
Transfected_only_without_outliers_1A <- Identify_and_remove_outliers(Reduced_transfected_1A_w_controls, 
                                                                     USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_1A_ap1_lowerbin <- Transfected_only_without_outliers_1A$ap1_lowerbin
Transfected_for_class_1A_ap1_upperbin <- Transfected_only_without_outliers_1A$ap1_upperbin
#Approach 2
Transfected_for_class_1A_ap2_lowerbin <- Transfected_only_without_outliers_1A$ap2_lowerbin
Transfected_for_class_1A_ap2_upperbin <- Transfected_only_without_outliers_1A$ap2_upperbin
#Approach 3
Transfected_for_class_1A_ap3_lowerbin <- Transfected_only_without_outliers_1A$ap3_lowerbin
Transfected_for_class_1A_ap3_upperbin <- Transfected_only_without_outliers_1A$ap3_upperbin


## Set 1A alternative comparisons
Transfected_only_without_outliers_1A_altcomps <- Identify_and_remove_outliers(Reduced_transfected_1Aagg_fused_ALTERNATIVE_COMPARISONS, 
                                                                     USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_1A_altcomps_ap1_lowerbin <- Transfected_only_without_outliers_1A_altcomps$ap1_lowerbin
Transfected_for_class_1A_altcomps_ap1_upperbin <- Transfected_only_without_outliers_1A_altcomps$ap1_upperbin
#Approach 2
Transfected_for_class_1A_altcomps_ap2_lowerbin <- Transfected_only_without_outliers_1A_altcomps$ap2_lowerbin
Transfected_for_class_1A_altcomps_ap2_upperbin <- Transfected_only_without_outliers_1A_altcomps$ap2_upperbin
#Approach 3
Transfected_for_class_1A_altcomps_ap3_lowerbin <- Transfected_only_without_outliers_1A_altcomps$ap3_lowerbin
Transfected_for_class_1A_altcomps_ap3_upperbin <- Transfected_only_without_outliers_1A_altcomps$ap3_upperbin


## Set B
Transfected_only_without_outliers_B <- Identify_and_remove_outliers(Reduced_transfected_B_w_controls, 
                                                                    USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_B_ap1_lowerbin <- Transfected_only_without_outliers_B$ap1_lowerbin
Transfected_for_class_B_ap1_upperbin <- Transfected_only_without_outliers_B$ap1_upperbin
#Approach 2
Transfected_for_class_B_ap2_lowerbin <- Transfected_only_without_outliers_B$ap2_lowerbin
Transfected_for_class_B_ap2_upperbin <- Transfected_only_without_outliers_B$ap2_upperbin
#Approach 3
Transfected_for_class_B_ap3_lowerbin <- Transfected_only_without_outliers_B$ap3_lowerbin
Transfected_for_class_B_ap3_upperbin <- Transfected_only_without_outliers_B$ap3_upperbin


## Set C
Transfected_only_without_outliers_C <- Identify_and_remove_outliers(Reduced_transfected_C_w_controls, 
                                                                    USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_C_ap1_lowerbin <- Transfected_only_without_outliers_C$ap1_lowerbin
Transfected_for_class_C_ap1_upperbin <- Transfected_only_without_outliers_C$ap1_upperbin
#Approach 2
Transfected_for_class_C_ap2_lowerbin <- Transfected_only_without_outliers_C$ap2_lowerbin
Transfected_for_class_C_ap2_upperbin <- Transfected_only_without_outliers_C$ap2_upperbin
#Approach 3
Transfected_for_class_C_ap3_lowerbin <- Transfected_only_without_outliers_C$ap3_lowerbin
Transfected_for_class_C_ap3_upperbin <- Transfected_only_without_outliers_C$ap3_upperbin


## Set D
Transfected_only_without_outliers_D <- Identify_and_remove_outliers(Reduced_transfected_D_w_controls, 
                                                                    USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_D_ap1_lowerbin <- Transfected_only_without_outliers_D$ap1_lowerbin
Transfected_for_class_D_ap1_upperbin <- Transfected_only_without_outliers_D$ap1_upperbin
#Approach 2
Transfected_for_class_D_ap2_lowerbin <- Transfected_only_without_outliers_D$ap2_lowerbin
Transfected_for_class_D_ap2_upperbin <- Transfected_only_without_outliers_D$ap2_upperbin
#Approach 3
Transfected_for_class_D_ap3_lowerbin <- Transfected_only_without_outliers_D$ap3_lowerbin
Transfected_for_class_D_ap3_upperbin <- Transfected_only_without_outliers_D$ap3_upperbin


## Set D alternative comparisons
Transfected_only_without_outliers_D_altcomps <- Identify_and_remove_outliers(Reduced_transfected_D_alt_comps, 
                                                                    USE_LOWER_MAX_UPPER_HINGE = TRUE)
#Approach 1
Transfected_for_class_D_altcomps_ap1_lowerbin <- Transfected_only_without_outliers_D_altcomps$ap1_lowerbin %>% 
  dplyr::select(-all_of(c('Intensity_MinIntensity_MCHERRY_BACKCORR', 'RadialDistribution_ZernikePhase_MCHERRY_BACKCORR_2_0')))
Transfected_for_class_D_altcomps_ap1_upperbin <- Transfected_only_without_outliers_D_altcomps$ap1_upperbin%>% 
  dplyr::select(-all_of(c('Intensity_MinIntensity_MCHERRY_BACKCORR', 'RadialDistribution_ZernikePhase_MCHERRY_BACKCORR_2_0')))
#Approach 2
Transfected_for_class_D_altcomps_ap2_lowerbin <- Transfected_only_without_outliers_D_altcomps$ap2_lowerbin%>% 
  dplyr::select(-all_of(c('Intensity_MinIntensity_MCHERRY_BACKCORR', 'RadialDistribution_ZernikePhase_MCHERRY_BACKCORR_2_0')))
Transfected_for_class_D_altcomps_ap2_upperbin <- Transfected_only_without_outliers_D_altcomps$ap2_upperbin%>% 
  dplyr::select(-all_of(c('Intensity_MinIntensity_MCHERRY_BACKCORR', 'RadialDistribution_ZernikePhase_MCHERRY_BACKCORR_2_0')))
#Approach 3
Transfected_for_class_D_altcomps_ap3_lowerbin <- Transfected_only_without_outliers_D_altcomps$ap3_lowerbin%>% 
  dplyr::select(-all_of(c('Intensity_MinIntensity_MCHERRY_BACKCORR', 'RadialDistribution_ZernikePhase_MCHERRY_BACKCORR_2_0')))
Transfected_for_class_D_altcomps_ap3_upperbin <- Transfected_only_without_outliers_D_altcomps$ap3_upperbin%>% 
  dplyr::select(-all_of(c('Intensity_MinIntensity_MCHERRY_BACKCORR', 'RadialDistribution_ZernikePhase_MCHERRY_BACKCORR_2_0')))




#

#NegativeTransfected cells for each isoform pair and controls
#Set 1
Negtrans_for_class_1 <- Identify_and_remove_outliers(Reduced_negtrans_1_w_controls, 
                                                     USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                     ONLY_NEGTRANS = TRUE,
                                                     REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                     TOP_x_PERCENTILE = 0.1,
                                                     REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                     TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_1_ap1_lowerbin,
                                                     relative_intensity_threshold_for_binning = 0.2,
                                                     relative_intensity_percentile_bins = 0.5)$full_dataset

#Set A
Negtrans_for_class_A <- Identify_and_remove_outliers(Reduced_negtrans_A_w_controls, 
                                                                       USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                                       ONLY_NEGTRANS = TRUE,
                                                                       REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                                       TOP_x_PERCENTILE = 0.1,
                                                                       REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                                       TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_A_ap1_lowerbin,
                                                                       relative_intensity_threshold_for_binning = 0.2,
                                                                       relative_intensity_percentile_bins = 0.5)$full_dataset

#Set 1A
Negtrans_for_class_1A <- Identify_and_remove_outliers(Reduced_negtrans_1A_w_controls, 
                                                                        USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                                        ONLY_NEGTRANS = TRUE,
                                                                        REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                                        TOP_x_PERCENTILE = 0.1,
                                                                        REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                                        TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_1A_ap1_lowerbin,
                                                                        relative_intensity_threshold_for_binning = 0.2,
                                                                        relative_intensity_percentile_bins = 0.5)$full_dataset



# Set 1A alternative comparisons
Negtrans_for_class_1A_altcomp <- Identify_and_remove_outliers(Reduced_negtrans_1Aagg_fused_ALTERNATIVE_COMPARISONS, 
                                                      USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                      ONLY_NEGTRANS = TRUE,
                                                      REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                      TOP_x_PERCENTILE = 0.1,
                                                      REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                      TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_1A_altcomps_ap1_lowerbin,
                                                      relative_intensity_threshold_for_binning = 0.2,
                                                      relative_intensity_percentile_bins = 0.5)$full_dataset

#Set B
Negtrans_for_class_B <- Identify_and_remove_outliers(Reduced_negtrans_B_w_controls, 
                                                                       USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                                       ONLY_NEGTRANS = TRUE,
                                                                       REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                                       TOP_x_PERCENTILE = 0.1,
                                                                       REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                                       TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_B_ap1_lowerbin,
                                                                       relative_intensity_threshold_for_binning = 0.2,
                                                                       relative_intensity_percentile_bins = 0.5)$full_dataset



#Set C
Negtrans_for_class_C <- Identify_and_remove_outliers(Reduced_negtrans_C_w_controls, 
                                                     USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                     ONLY_NEGTRANS = TRUE,
                                                     REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                     TOP_x_PERCENTILE = 0.1,
                                                     REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                     TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_C_ap1_lowerbin,
                                                     relative_intensity_threshold_for_binning = 0.2,
                                                     relative_intensity_percentile_bins = 0.5)$full_dataset




#Set D
Negtrans_for_class_D <- Identify_and_remove_outliers(Reduced_negtrans_D_w_controls, 
                                                     USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                     ONLY_NEGTRANS = TRUE,
                                                     REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                     TOP_x_PERCENTILE = 0.1,
                                                     REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                     TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_D_ap1_lowerbin,
                                                     relative_intensity_threshold_for_binning = 0.2,
                                                     relative_intensity_percentile_bins = 0.5)$full_dataset

#Set D alt comp
#Set D
Negtrans_for_class_D_altcomp <- Identify_and_remove_outliers(Reduced_negtrans_D_alt_comps, 
                                                     USE_LOWER_MAX_UPPER_HINGE = TRUE,
                                                     ONLY_NEGTRANS = TRUE,
                                                     REMOVE_TOP_x_PERCENTILE_NEGATIVE = TRUE,
                                                     TOP_x_PERCENTILE = 0.1,
                                                     REMOVE_SAME_COLS_AS_TRANSFECTED = TRUE,
                                                     TRANSFECTED_DATASET_FOR_COLS_TO_KEEP = Transfected_for_class_D_altcomps_ap1_lowerbin,
                                                     relative_intensity_threshold_for_binning = 0.2,
                                                     relative_intensity_percentile_bins = 0.5)$full_dataset

#####

#CHECK THAT THE SAME COLUMNS ARE USED FOR CLASSIFICATION BETWEEN NON-TRANSFECTED AND POSITIVE TRANSFECTANTS
#####
#Set 1
all(colnames(Negtrans_for_class_1) == colnames(Transfected_for_class_1_ap1_lowerbin))
#Set A
all(colnames(Negtrans_for_class_A) == colnames(Transfected_for_class_A_ap1_lowerbin))
#Set 1A
all(colnames(Negtrans_for_class_1A) == colnames(Transfected_for_class_1A_ap1_lowerbin))
#Set 1A altcomp
all(colnames(Negtrans_for_class_1A_altcomp) == colnames(Transfected_for_class_1A_altcomps_ap1_lowerbin))
#Set B
all(colnames(Negtrans_for_class_B) == colnames(Transfected_for_class_B_ap1_lowerbin))
#Set C
all(colnames(Negtrans_for_class_C) == colnames(Transfected_for_class_C_ap1_lowerbin))
#Set D
all(colnames(Negtrans_for_class_D) == colnames(Transfected_only_without_outliers_D))
#Set D altcomp
all(colnames(Negtrans_for_class_D_altcomp) == colnames(Transfected_for_class_D_altcomps_ap1_lowerbin))



# Function to sample dataset X based on proportions in dataset Y
# While maintaining equal numbers of each ISOFORM value for each replicate
# If insufficient observations in any replicate/isoform combination, 
# downsamples all replicates proportionally
#
# X: Main dataset to sample from
# Y: Dataset containing the counts for each Experiment_Rep
# seed: Optional random seed for reproducibility
proportional_sampling <- function(X, Y, seed = NULL) {
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Error handling - check if datasets exist
  if (missing(X) || is.null(X) || nrow(X) == 0) {
    warning("Dataset X is missing, NULL, or empty")
    return(NULL)
  }
  
  if (missing(Y) || is.null(Y) || nrow(Y) == 0) {
    warning("Dataset Y (counts) is missing, NULL, or empty")
    return(NULL)
  }
  
  # Check if required columns exist in both datasets
  if (!"Experiment_Rep" %in% colnames(X)) {
    warning("Column 'Experiment_Rep' not found in dataset X")
    return(NULL)
  }
  
  if (!"ISOFORM" %in% colnames(X)) {
    warning("Column 'ISOFORM' not found in dataset X")
    return(NULL)
  }
  
  if (!"Experiment_Rep" %in% colnames(Y)) {
    warning("Column 'Experiment_Rep' not found in dataset Y")
    return(NULL)
  }
  
  if (!"n" %in% colnames(Y)) {
    warning("Column 'n' not found in dataset Y")
    return(NULL)
  }
  
  # Check if ISOFORM column contains exactly 2 unique values
  isoform_values <- unique(X$ISOFORM)
  if (length(isoform_values) != 2) {
    warning("Dataset X must contain exactly 2 unique values in the ISOFORM column, found: ", 
            paste(isoform_values, collapse = ", "))
    return(NULL)
  }
  
  # Get the unique replicates in X
  replicates_in_X <- unique(X$Experiment_Rep)
  
  # Filter Y to only include replicates that exist in X
  Y_filtered <- Y[Y$Experiment_Rep %in% replicates_in_X, ]
  
  if (nrow(Y_filtered) == 0) {
    warning("No matching Experiment_Rep values found between datasets")
    return(NULL)
  }
  
  # First pass: Check all replicate and isoform combinations to find the limiting one
  rep_iso_counts <- data.frame()
  
  # Calculate desired counts and available counts for each rep/isoform combination
  for (rep in Y_filtered$Experiment_Rep) {
    # Extract the target count for this replicate
    target_count <- Y_filtered$n[Y_filtered$Experiment_Rep == rep]
    
    # Since we need equal numbers of each ISOFORM, ensure target_count is even
    if (target_count %% 2 != 0) {
      target_count <- target_count + 1
    }
    
    # Number to sample from each ISOFORM type
    per_isoform_count <- target_count / 2
    
    # Check available observations for each isoform
    for (iso in isoform_values) {
      available_count <- nrow(X[X$Experiment_Rep == rep & X$ISOFORM == iso, ])
      rep_iso_counts <- rbind(rep_iso_counts, data.frame(
        Experiment_Rep = rep,
        ISOFORM = iso,
        desired_count = per_isoform_count,
        available_count = available_count,
        ratio = available_count / per_isoform_count
      ))
    }
  }
  
  # Find the limiting ratio (the smallest ratio of available/desired)
  min_ratio <- min(rep_iso_counts$ratio)
  
  # If min_ratio < 1, we need to downsample everything proportionally
  if (min_ratio < 1) {
    warning(paste("Insufficient observations detected. Downsampling all replicates by factor:", 
                  round(min_ratio, 3)))
    
    # Adjust all desired counts by the min_ratio
    rep_iso_counts$adjusted_count <- floor(rep_iso_counts$desired_count * min_ratio)
    
    # Ensure we have at least 1 observation per rep/isoform combination
    rep_iso_counts$adjusted_count <- pmax(rep_iso_counts$adjusted_count, 1)
  } else {
    # If we have enough observations, use the original desired counts
    rep_iso_counts$adjusted_count <- rep_iso_counts$desired_count
  }
  
  # Second pass: Sample according to adjusted counts
  result <- NULL
  
  for (i in 1:nrow(rep_iso_counts)) {
    rep <- rep_iso_counts$Experiment_Rep[i]
    iso <- rep_iso_counts$ISOFORM[i]
    adjusted_count <- rep_iso_counts$adjusted_count[i]
    
    # Get data for this replicate and isoform
    rep_iso_data <- X[X$Experiment_Rep == rep & X$ISOFORM == iso, ]
    
    # Sample without replacement (we know we have enough now)
    sampled_data <- rep_iso_data[sample(1:nrow(rep_iso_data), adjusted_count, replace = FALSE), ]
    
    # Append to result
    result <- rbind(result, sampled_data)
  }
  
  # Calculate the actual proportions achieved to inform the user
  result_summary <- as.data.frame(table(result$Experiment_Rep))
  colnames(result_summary) <- c("Experiment_Rep", "actual_count")
  
  Y_summary <- as.data.frame(tapply(Y_filtered$n, Y_filtered$Experiment_Rep, sum))
  Y_summary$Experiment_Rep <- rownames(Y_summary)
  colnames(Y_summary)[1] <- "desired_count"
  
  proportions_summary <- merge(Y_summary, result_summary, by = "Experiment_Rep")
  proportions_summary$desired_prop <- proportions_summary$desired_count / sum(proportions_summary$desired_count)
  proportions_summary$actual_prop <- proportions_summary$actual_count / sum(proportions_summary$actual_count)
  proportions_summary$difference <- proportions_summary$actual_prop - proportions_summary$desired_prop
  
  # Print summary
  message("Sampling summary:")
  message("Min ratio of available/desired observations: ", round(min_ratio, 3))
  message("Total observations sampled: ", nrow(result))
  message("Replicate proportion differences (actual - desired):")
  for (i in 1:nrow(proportions_summary)) {
    message(paste0("  ", proportions_summary$Experiment_Rep[i], ": ", 
                   sprintf("%.3f", proportions_summary$difference[i])))
  }
  
  # Verify that we have equal numbers of each isoform within each replicate
  iso_balance <- as.data.frame(table(result$Experiment_Rep, result$ISOFORM))
  colnames(iso_balance) <- c("Experiment_Rep", "ISOFORM", "Count")
  iso_balance_wide <- reshape(iso_balance, 
                              idvar = "Experiment_Rep",
                              timevar = "ISOFORM", 
                              direction = "wide")
  
  if (ncol(iso_balance_wide) > 2) {  # Check if reshape worked correctly
    imbalance <- any(abs(iso_balance_wide[,2] - iso_balance_wide[,3]) > 0)
    if (imbalance) {
      warning("Isoform classes are not perfectly balanced within some replicates")
    } else {
      message("Isoform classes are perfectly balanced within each replicate")
    }
  }
  
  return(result)
}




### Classifier ### 

#
isoform_classify <- function(DATASET,
                             PAIR,
                             FEATURE_START = 14,
                             model_type = 'nnet',
                             downsample_for_equal_classes,
                             sample_controls_to_match_pair,
                             include_positive_control,
                             MCHERRY_thresholding_with_max_mean_when_downsampling,
                             reduce_model_output,
                             print_training_progress = FALSE,
                             include_negative_control){
  ### Downsampling isoforms
  if(downsample_for_equal_classes == TRUE){
    ISO_DF <- DATASET %>% 
      dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == PAIR)
    iso_set_pair <- ISO_DF %>% 
      dplyr::pull(SET_PAIR) %>% 
      unique()
    ##### 
    #n cells per iso per rep
    iso1 <- ISO_DF %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
      summarize(n_cells = n(),
                'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
    #Get numbers of cells for each isoform in a given pair with least no. of cells <- this is the n_cells to sample for a given gene in a given rep
    iso2 <- iso1 %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep) %>%
      summarize('n_cells_to_sample' = min(n_cells),
                'CHERRY_threshold' = min(max_MCHERRY_intensity))
    #Join tibble with n_cells_to_sammple to original dataset to get sample no. fused w. SET_PAIR_REP_ISOFORM
    iso3 <- left_join(x = ISO_DF, y = iso2, by = c('SET_PAIR','Experiment_Rep')) %>% 
      dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
      unique()
    #This codes sammples the n_cells_to_sample for each isoform pair in a given rep (note, slice_sample() only takes a single value input and as we want to sample different numbers per group, this is what we need to do)
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
      iso_a <- try(
        ISO_DF %>%
          left_join(x = ISO_DF, y = iso3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
          dplyr::select(!c(n_cells_to_sample)) %>%
          dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n),silent = TRUE)
      if(inherits(iso_a, "try-error")){
        iso4 <- ISO_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
        iso_cherry_match <- FALSE
      }
      else{
        iso4 <- iso_a
        iso_cherry_match <- TRUE
      }
    }
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
      iso4 <- ISO_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR_REP_ISOFORM) %>% 
        nest() %>%            
        ungroup() %>% 
        mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
        mutate(samp = map2(data, n, sample_n)) %>% 
        select(-data) %>%
        unnest(samp) %>% 
        dplyr::select(!n)
    }
    #Count
    iso5 <- iso4 %>% 
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      count()
    #
    ISO_DATASET_TO_TRIM <- iso4 #
    
    #
    n_postrans_for_balanced_control <- iso2 %>% 
      dplyr::select(SET_PAIR, Experiment_Rep, n_cells_to_sample) %>% 
      dplyr::mutate(n = n_cells_to_sample * 2) %>% 
      dplyr::select(-n_cells_to_sample)
  }
  #####
  ### If include positive_control is true, this generates the needed objects
  #####
  if(include_positive_control == TRUE){
    #
    POSITIVE_DF <- DATASET %>%
      ungroup() %>% 
      dplyr::filter(c(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A' | word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_MCHERRY')) %>% 
      dplyr::mutate(ISOFORM = case_when(
        ISOFORM %in% c('WT','DEx') ~ 'RAB5A',
        ISOFORM == 'MCHERRY' ~ 'MCHERRY'),
        SET_PAIR = paste0(Experiment_Set,'_PAIR_POSITIVE'),
        SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
    pos_set_pair <- POSITIVE_DF %>% 
      pull(SET_PAIR) %>% 
      unique()
    #
    POS1 <- POSITIVE_DF  %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
      summarize(n_cells = n(),
                'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
    #
    POS2 <- POS1 %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep) %>%
      summarize('n_cells_to_sample' = min(n_cells),
                'CHERRY_threshold' = min(max_MCHERRY_intensity))
    #
    POS3 <- left_join(x = POSITIVE_DF, y = POS2, by = c('SET_PAIR','Experiment_Rep')) %>% 
      dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
      unique()
    #
    
    #
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
      pos_a <- try(
        POSITIVE_DF %>%
          left_join(x = POSITIVE_DF, y = POS3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
          dplyr::select(!c(n_cells_to_sample)) %>%
          dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% POS3$SET_PAIR_REP_ISOFORM, POS3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n),
        silent = TRUE)
      # }
      if(inherits(pos_a, "try-error") | iso_cherry_match == FALSE){
        POS4 <- POSITIVE_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% POS3$SET_PAIR_REP_ISOFORM, POS3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
        pos_cherry_match <- FALSE
      }
      if(inherits(pos_a, 'tbl') & iso_cherry_match == TRUE){
        POS4 <- pos_a
        pos_cherry_match <- TRUE
      }
    }
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
      POS4 <- POSITIVE_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR_REP_ISOFORM) %>% 
        nest() %>%            
        ungroup() %>% 
        mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% POS3$SET_PAIR_REP_ISOFORM, POS3$n_cells_to_sample, 0)) %>%
        mutate(samp = map2(data, n, sample_n)) %>% 
        select(-data) %>%
        unnest(samp) %>% 
        dplyr::select(!n)
    }
    #
    POS5 <- POS4 %>% 
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      count()
    
  }
  #####
  ### Downsampling RAB5A to match isoform (if downsampling of isoform has happened or not)
  #####
  if(include_negative_control == TRUE){
    if(downsample_for_equal_classes == TRUE && sample_controls_to_match_pair == TRUE){
      #Isoform count
      n_iso <- iso5 %>% 
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == PAIR) %>% 
        dplyr::pull(n)
      #Negative
      NEG_DF <- DATASET %>%
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A')
      
      neg_set_pair <- NEG_DF %>% 
        dplyr::pull(SET_PAIR) %>% 
        unique()
      #n cells per iso per rep
      NEG1 <- NEG_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
        summarize(n_cells = n(),
                  'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
      #Get numbers of cells for each isoform in a given pair with least no. of cells <- this is the n_cells to sample for a given gene in a given rep
      NEG2 <- NEG1 %>% 
        ungroup() %>% 
        group_by(SET_PAIR, Experiment_Rep) %>%
        summarize('n_cells_to_sample' = min(n_cells),
                  'CHERRY_threshold' = min(max_MCHERRY_intensity))
      #Join tibble with n_cells_to_sammple to original dataset to get sample no. fused w. SET_PAIR_REP_ISOFORM
      NEG3 <- left_join(x = NEG_DF, y = NEG2, by = c('SET_PAIR','Experiment_Rep')) %>% 
        dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
        unique()
      #This codes sammples the n_cells_to_sample for each isoform pair in a given rep (note, slice_sample() only takes a single value input and as we want to sample different numbers per group, this is what we need to do)
      if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
        neg_a <- try(
          NEG_DF %>%
            left_join(x = NEG_DF, y = NEG3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
            dplyr::select(!c(n_cells_to_sample)) %>%
            dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
            ungroup() %>% 
            group_by(SET_PAIR_REP_ISOFORM) %>% 
            nest() %>%            
            ungroup() %>% 
            mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% NEG3$SET_PAIR_REP_ISOFORM, NEG3$n_cells_to_sample, NA)) %>%
            mutate(samp = map2(data, n, sample_n)) %>% 
            select(-data) %>%
            unnest(samp) %>% 
            dplyr::select(!n),silent = TRUE)
        
        if(inherits(neg_a, "try-error") | iso_cherry_match == FALSE){
          NEG4 <- NEG_DF %>% 
            ungroup() %>% 
            group_by(SET_PAIR_REP_ISOFORM) %>% 
            nest() %>%            
            ungroup() %>% 
            mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% NEG3$SET_PAIR_REP_ISOFORM, NEG3$n_cells_to_sample, NA)) %>%
            mutate(samp = map2(data, n, sample_n)) %>% 
            select(-data) %>%
            unnest(samp) %>% 
            dplyr::select(!n)
          neg_cherry_match <- FALSE
        }
        if(inherits(neg_a, 'tbl') & iso_cherry_match == TRUE){
          NEG4 <- neg_a
          neg_cherry_match <- TRUE
        }
      }
      if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
        NEG4 <- NEG_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% NEG3$SET_PAIR_REP_ISOFORM, NEG3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
      }
      #Negative control Count
      n_negative <- NEG4 %>% 
        ungroup() %>% 
        group_by(SET_PAIR) %>% 
        count() %>% 
        dplyr::pull(n)
      
      
      
      #
      NEGATIVE_DATASET_TO_TRIM <- proportional_sampling(X =  NEG4,
                                                        Y = n_postrans_for_balanced_control,
                                                        seed = NULL)
    }
      if(include_positive_control == TRUE){
        #
        POSITIVE_DATASET_TO_TRIM <- proportional_sampling(X =  POS4,
                                                          Y = n_postrans_for_balanced_control,
                                                          seed = NULL)
      }
  }
  ###
  if(downsample_for_equal_classes == FALSE && sample_controls_to_match_pair == TRUE){
    n_iso <- DATASET %>% 
      dplyr::mutate(GENE_PAIR = word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))) %>% 
      dplyr::filter(GENE_PAIR == PAIR) %>%
      dim()
    n_iso <- n_iso[1]
    
    n_negative <- DATASET %>% 
      dplyr::mutate(GENE_PAIR = word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))) %>% 
      dplyr::filter(GENE_PAIR == 'PAIR_RAB5A') %>%
      dim()
    n_negative <- n_negative[1]
    if(n_negative > n_iso){
      NEGATIVE_DATASET_TO_TRIM <- DATASET %>%
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A') %>% 
        slice_sample(n = n_iso)
    }
    if(n_negative <= n_iso){
      NEGATIVE_DATASET_TO_TRIM <- DATASET %>% 
        ungroup() %>%
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A')
    }
    if(include_positive_control == TRUE){
      n_positive <- DATASET %>%
        ungroup() %>% 
        dplyr::filter(word(SET_PAIR, star = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A' | word(SET_PAIR, star = 3, end = 4, sep = fixed('_')) == 'PAIR_MCHERRY') %>% 
        dplyr::mutate(ISOFORM = case_when(
          ISOFORM %in% c('WT','DEx') ~ 'RAB5A',
          ISOFORM == 'MCHERRY' ~ 'MCHERRY'),
          SET_PAIR = paste0(Experiment_Set,'_PAIR_POSITIVE'),
          SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM)) %>% 
        dim()
      n_positive <- n_positive[1]
      
      if(n_positive > n_iso){
        POSITIVE_DATASET_TO_TRIM <- POSITIVE_DF %>% 
          slice_sample(n = n_iso)
      }
      if(n_positive <= n_iso){
        POSITIVE_DATASET_TO_TRIM <- POSITIVE_DF
      }
    }
  }
  #####
  ### Downsampling positive control (currently using RAB5A cells from both wells)
  #####
  
  #####
  ### Creating variable with correct name if downsampling was not done
  if(downsample_for_equal_classes == FALSE){
    ISO_DATASET_TO_TRIM <- DATASET
  }
  ###
  ISOFORM_DATASET <- ISO_DATASET_TO_TRIM %>% 
    dplyr::mutate(GENE_PAIR = word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))) %>% 
    dplyr::filter(GENE_PAIR == PAIR) %>%
    dplyr::select(ISOFORM, all_of(FEATURE_START:length(colnames(DATASET)))) %>%
    dplyr::select(!c(starts_with('Intensity') & ends_with('MCHERRY_BACKCORR') & !contains('MeanIntensity_'))) #Removes intensity data (although leaves intensities normalised to mean intensities)
  #
  if(include_negative_control == TRUE){
    NEGATIVE_DATASET <- NEGATIVE_DATASET_TO_TRIM %>% 
      dplyr::select(ISOFORM, all_of(FEATURE_START:length(colnames(DATASET)))) %>%
      dplyr::select(!c(starts_with('Intensity') & ends_with('MCHERRY_BACKCORR') & !contains('MeanIntensity_'))) #Removes intensity data (although leaves intensities normalised to mean intensities)
  }
  # 
  if(include_positive_control == TRUE){
    POSITIVE_DATASET <- POSITIVE_DATASET_TO_TRIM %>% 
      dplyr::select(ISOFORM, all_of(FEATURE_START:length(colnames(DATASET)))) %>%
      dplyr::select(!c(starts_with('Intensity') & ends_with('MCHERRY_BACKCORR') & !contains('MeanIntensity_'))) #Removes intensity data (although leaves intensities normalised to mean intensities)
  }
  #
  if(include_positive_control == FALSE & include_negative_control == FALSE){
    dataset_list <- list('ISOFORM_DF' = ISOFORM_DATASET)
  }
  #
  if(include_positive_control == FALSE & include_negative_control == TRUE){
    dataset_list <- list('ISOFORM_DF' = ISOFORM_DATASET, 'NEGATIVE_DF' = NEGATIVE_DATASET)
  }
  #
  if(include_positive_control == TRUE & include_negative_control == TRUE){
    dataset_list <- list('ISOFORM_DF' = ISOFORM_DATASET, 'NEGATIVE_DF' = NEGATIVE_DATASET, 'POSITIVE_DF' = POSITIVE_DATASET)
  }
  #return(dataset_list)
  
  
  ## Remove when done implementing
  #Caret models
  #####
  #Control object for caret models.   
  if(length(unique(ISOFORM_DATASET$ISOFORM)) == 2){
    myControl <- trainControl(
      method = "cv", #cv = cross validation
      number = 10, #Number of folds
      summaryFunction = twoClassSummary,
      classProbs = TRUE, # <- Super important!
      verboseIter = print_training_progress,
      preProcOptions = list(thresh = 0.9, freqCut = 70/30, uniqueCut = 25), #Note, this is the way to modify how the arguments in the proprocess function works when used within the train workflow
    )
  }
  if(length(unique(ISOFORM_DATASET$ISOFORM)) > 2){
    myControl <- trainControl(
      method = "cv", #cv = cross validation
      number = 10, #Number of folds
      summaryFunction = defaultSummary,
      classProbs = TRUE, # <- Super important!
      verboseIter = print_training_progress,
      preProcOptions = list(thresh = 0.9, freqCut = 70/30, uniqueCut = 25), #Note, this is the way to modify how the arguments in the proprocess function works when used within the train workflow
    )
  }
  
  if(model_type == 'glmnet'){
    #Training glmnet model
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      preProcess = c('zv','nzv','center','scale'),
      trControl = myControl
    )
    #  negative_model <- train(
    #    ISOFORM ~ .,
    #    data = NEGATIVE_DATASET,
    #    method = model_type,
    #    preProcess = c('zv','nzv','center','scale'),
    #    trControl = myControl
    #  )
    #  if(include_positive_control == TRUE){
    #    positive_model <- train(
    #      ISOFORM ~ .,
    #      data = POSITIVE_DATASET,
    #      method = model_type,
    #      preProcess = c('zv','nzv','center','scale'),
    #      trControl = myControl
    #    )
    # }
  }
  #Svm classifiers - can be extended if you want other kernels etc.
  if(model_type %in% c('svmLinear2')){
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      preProcess = c('zv','center','scale','pca'),
      trControl = myControl
    )
  }
  #Ranger (a random-forest type)
  if(model_type %in% c('ranger','rf')){
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      preProcess = c('zv','center','scale'),
      trControl = myControl
    )
    negative_model <- train(
      ISOFORM ~ .,
      data = NEGATIVE_DATASET,
      method = model_type,
      preProcess = c('zv','center','scale'),
      trControl = myControl
    )
    if(include_positive_control == TRUE){
      positive_model <- train(
        ISOFORM ~ .,
        data = POSITIVE_DATASET,
        method = model_type,
        preProcess = c('zv','center','scale'),
        trControl = myControl
      )
    }
  }
  #Neural network (feedforward neural network, FNN)
  if(model_type == 'nnet'){
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      MaxNWts = 1060,
      preProcess = c('zv','center','scale'),
      trControl = myControl
    )
    if(include_negative_control == TRUE){
      negative_model <- train(
        ISOFORM ~ .,
        data = NEGATIVE_DATASET,
        method = model_type,
        MaxNWts = 1060,
        preProcess = c('zv','center','scale'),
        trControl = myControl
      )
    }
    if(include_positive_control == TRUE){
      positive_model <- train(
        ISOFORM ~ .,
        data = POSITIVE_DATASET,
        method = model_type,
        MaxNWts = 1060,
        preProcess = c('zv','center','scale'),
        trControl = myControl
      )
    }
  }
  #####
  if(reduce_model_output == TRUE){
    reduced_model_lists <- function(LIST){
      #
      variable_importance <- varImp(LIST)$importance
      #
      LIST$trainingData <- LIST$trainingData %>% 
        dplyr::select(.outcome, Intensity_MeanIntensity_MCHERRY_BACKCORR)
      #
      reduced_list <- list('results' = LIST$results, 'training_data_intensity' = LIST$trainingData, 'variable_importance' = variable_importance)
      return(reduced_list)
    }
    #
    isoform_model <- reduced_model_lists(isoform_model)
    isoform_model <- append(isoform_model, values = c(list('iso_cherry_match' = iso_cherry_match),list('SET_PAIR' = iso_set_pair)))
    
    if(include_negative_control == TRUE){
      negative_model <- reduced_model_lists(negative_model)
      negative_model <- append(negative_model, values = c(list('neg_cherry_match' = neg_cherry_match),list('SET_PAIR' = neg_set_pair)))
    }
    #
    if(include_positive_control == TRUE){
      positive_model <- reduced_model_lists(positive_model)
      positive_model <- append(positive_model, values = c(list('pos_cherry_match' = pos_cherry_match),list('SET_PAIR' = pos_set_pair)))
    }
  }
  if(include_positive_control == FALSE & include_negative_control == FALSE){
    model_list <- list('ISOFORM_MODEL' = isoform_model)
  }
  if(include_positive_control == FALSE & include_negative_control == TRUE){
    model_list <- list('ISOFORM_MODEL' = isoform_model, 'NEGATIVE_MODEL' = negative_model)
  }
  if(include_positive_control == TRUE & include_negative_control == TRUE){
    model_list <- list('ISOFORM_MODEL' = isoform_model, 'NEGATIVE_MODEL' = negative_model, 'POSITIVE_MODEL' = positive_model)
  }
  return(model_list)
}

#
isoform_classify_for_negative_cells <- function(DATASET,
                                                PAIR,
                                                FEATURE_START = 14,
                                                model_type = 'nnet',
                                                downsample_for_equal_classes,
                                                sample_controls_to_match_pair,
                                                include_positive_control,
                                                MCHERRY_thresholding_with_max_mean_when_downsampling,
                                                reduce_model_output,
                                                print_training_progress = FALSE,
                                                POSITIVE_DATASET,
                                                include_negative_control){
  ### Downsampling isoforms + generating count for POSITIVE TRANSFECTANTS
  #We will use this to reduce sample size of NEGATIVE TRANSFECTANTS so they are matched
  if(downsample_for_equal_classes == TRUE){
    ISO_DF <- POSITIVE_DATASET %>% 
      dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == PAIR)
    iso_set_pair <- ISO_DF %>% 
      dplyr::pull(SET_PAIR) %>% 
      unique()
    ##### 
    #n cells per iso per rep
    iso1 <- ISO_DF %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
      summarize(n_cells = n(),
                'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
    #Get numbers of cells for each isoform in a given pair with least no. of cells <- this is the n_cells to sample for a given gene in a given rep
    iso2 <- iso1 %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep) %>%
      summarize('n_cells_to_sample' = min(n_cells),
                'CHERRY_threshold' = min(max_MCHERRY_intensity))
    #Join tibble with n_cells_to_sammple to original dataset to get sample no. fused w. SET_PAIR_REP_ISOFORM
    iso3 <- left_join(x = ISO_DF, y = iso2, by = c('SET_PAIR','Experiment_Rep')) %>% 
      dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
      unique()
    #This codes sammples the n_cells_to_sample for each isoform pair in a given rep (note, slice_sample() only takes a single value input and as we want to sample different numbers per group, this is what we need to do)
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
      iso_a <- try(
        ISO_DF %>%
          left_join(x = ISO_DF, y = iso3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
          dplyr::select(!c(n_cells_to_sample)) %>%
          dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n),silent = TRUE)
      if(inherits(iso_a, "try-error")){
        iso4 <- ISO_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
        iso_cherry_match <- FALSE
      }
      else{
        iso4 <- iso_a
        iso_cherry_match <- TRUE
      }
    }
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
      iso4 <- ISO_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR_REP_ISOFORM) %>% 
        nest() %>%            
        ungroup() %>% 
        mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
        mutate(samp = map2(data, n, sample_n)) %>% 
        select(-data) %>%
        unnest(samp) %>% 
        dplyr::select(!n)
      iso_cherry_match <- FALSE
    }
    #Count
    iso5 <- iso4 %>% 
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      count()
    #
    #ISO_DATASET_TO_TRIM <- iso4 #
    #Isoform count
    n_iso_POSITIVE <- iso5 %>% 
      dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == PAIR) %>% 
      dplyr::pull(n)
    
    #
    n_postrans_for_balanced_control <- iso2 %>% 
      dplyr::select(SET_PAIR, Experiment_Rep, n_cells_to_sample) %>% 
      dplyr::mutate(n = n_cells_to_sample * 2) %>% 
      dplyr::select(-n_cells_to_sample)
  }
  
  ### Downsampling isoforms
  if(downsample_for_equal_classes == TRUE){
    ISO_DF <- DATASET %>% 
      dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == PAIR)
    iso_set_pair <- ISO_DF %>% 
      dplyr::pull(SET_PAIR) %>% 
      unique()
    ##### 
    #n cells per iso per rep
    iso1 <- ISO_DF %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
      summarize(n_cells = n(),
                'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
    #Get numbers of cells for each isoform in a given pair with least no. of cells <- this is the n_cells to sample for a given gene in a given rep
    iso2 <- iso1 %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep) %>%
      summarize('n_cells_to_sample' = min(n_cells),
                'CHERRY_threshold' = min(max_MCHERRY_intensity))
    #Join tibble with n_cells_to_sammple to original dataset to get sample no. fused w. SET_PAIR_REP_ISOFORM
    iso3 <- left_join(x = ISO_DF, y = iso2, by = c('SET_PAIR','Experiment_Rep')) %>% 
      dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
      unique()
    #This codes sammples the n_cells_to_sample for each isoform pair in a given rep (note, slice_sample() only takes a single value input and as we want to sample different numbers per group, this is what we need to do)
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
      iso_a <- try(
        ISO_DF %>%
          left_join(x = ISO_DF, y = iso3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
          dplyr::select(!c(n_cells_to_sample)) %>%
          dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n),silent = TRUE)
      if(inherits(iso_a, "try-error")){
        iso4 <- ISO_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
        iso_cherry_match <- FALSE
      }
      else{
        iso4 <- iso_a
        iso_cherry_match <- TRUE
      }
    }
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
      iso4 <- ISO_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR_REP_ISOFORM) %>% 
        nest() %>%            
        ungroup() %>% 
        mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% iso3$SET_PAIR_REP_ISOFORM, iso3$n_cells_to_sample, NA)) %>%
        mutate(samp = map2(data, n, sample_n)) %>% 
        select(-data) %>%
        unnest(samp) %>% 
        dplyr::select(!n)
      iso_cherry_match <- FALSE
    }
    #Count
    iso5 <- iso4 %>% 
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      count()
    #
    ISO_DATASET_TO_TRIM <- iso4 #
    ISO_DATASET_TO_TRIM <- proportional_sampling(X =  iso4,
                                                 Y = n_postrans_for_balanced_control,
                                                 seed = NULL)
    
  }
  #####
  ### If include positive_control is true, this generates the needed objects
  #####
  if(include_positive_control == TRUE){
    #
    POSITIVE_DF <- DATASET %>%
      ungroup() %>% 
      dplyr::filter(c(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A' | word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_MCHERRY')) %>% 
      dplyr::mutate(ISOFORM = case_when(
        ISOFORM %in% c('WT','DEx') ~ 'RAB5A',
        ISOFORM == 'MCHERRY' ~ 'MCHERRY'),
        SET_PAIR = paste0(Experiment_Set,'_PAIR_POSITIVE'),
        SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM))
    pos_set_pair <- POSITIVE_DF %>% 
      pull(SET_PAIR) %>% 
      unique()
    #
    POS1 <- POSITIVE_DF  %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
      summarize(n_cells = n(),
                'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
    #
    POS2 <- POS1 %>% 
      ungroup() %>% 
      group_by(SET_PAIR, Experiment_Rep) %>%
      summarize('n_cells_to_sample' = min(n_cells),
                'CHERRY_threshold' = min(max_MCHERRY_intensity))
    #
    POS3 <- left_join(x = POSITIVE_DF, y = POS2, by = c('SET_PAIR','Experiment_Rep')) %>% 
      dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
      unique()
    #
    
    #
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
      pos_a <- try(
        POSITIVE_DF %>%
          left_join(x = POSITIVE_DF, y = POS3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
          dplyr::select(!c(n_cells_to_sample)) %>%
          dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% POS3$SET_PAIR_REP_ISOFORM, POS3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n),
        silent = TRUE)
      # }
      if(inherits(pos_a, "try-error") | iso_cherry_match == FALSE){
        POS4 <- POSITIVE_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% POS3$SET_PAIR_REP_ISOFORM, POS3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
        pos_cherry_match <- FALSE
      }
      if(inherits(pos_a, 'tbl') & iso_cherry_match == TRUE){
        POS4 <- pos_a
        pos_cherry_match <- TRUE
      }
    }
    if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
      POS4 <- POSITIVE_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR_REP_ISOFORM) %>% 
        nest() %>%            
        ungroup() %>% 
        mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% POS3$SET_PAIR_REP_ISOFORM, POS3$n_cells_to_sample, 0)) %>%
        mutate(samp = map2(data, n, sample_n)) %>% 
        select(-data) %>%
        unnest(samp) %>% 
        dplyr::select(!n)
      pos_cherry_match <- TRUE
    }
    #
    POS5 <- POS4 %>% 
      ungroup() %>% 
      group_by(SET_PAIR) %>% 
      count()
    
  }
  #####
  ### Downsampling RAB5A to match isoform (if downsampling of isoform has happened or not)
  #####
  if(include_negative_control == TRUE){
    if(downsample_for_equal_classes == TRUE && sample_controls_to_match_pair == TRUE){
      #Isoform count
      n_iso <- iso5 %>% 
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == PAIR) %>% 
        dplyr::pull(n)
      #Negative
      NEG_DF <- DATASET %>%
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A')
      
      neg_set_pair <- NEG_DF %>% 
        dplyr::pull(SET_PAIR) %>% 
        unique()
      #n cells per iso per rep
      NEG1 <- NEG_DF %>% 
        ungroup() %>% 
        group_by(SET_PAIR, Experiment_Rep, SET_PAIR_REP_ISOFORM) %>% 
        summarize(n_cells = n(),
                  'max_MCHERRY_intensity' = max(Intensity_MeanIntensity_MCHERRY_BACKCORR))
      #Get numbers of cells for each isoform in a given pair with least no. of cells <- this is the n_cells to sample for a given gene in a given rep
      NEG2 <- NEG1 %>% 
        ungroup() %>% 
        group_by(SET_PAIR, Experiment_Rep) %>%
        summarize('n_cells_to_sample' = min(n_cells),
                  'CHERRY_threshold' = min(max_MCHERRY_intensity))
      #Join tibble with n_cells_to_sammple to original dataset to get sample no. fused w. SET_PAIR_REP_ISOFORM
      NEG3 <- left_join(x = NEG_DF, y = NEG2, by = c('SET_PAIR','Experiment_Rep')) %>% 
        dplyr::select(SET_PAIR_REP_ISOFORM, n_cells_to_sample, CHERRY_threshold) %>% 
        unique()
      #This codes sammples the n_cells_to_sample for each isoform pair in a given rep (note, slice_sample() only takes a single value input and as we want to sample different numbers per group, this is what we need to do)
      if(MCHERRY_thresholding_with_max_mean_when_downsampling == TRUE){
        neg_a <- try(
          NEG_DF %>%
            left_join(x = NEG_DF, y = NEG3, by = c('SET_PAIR_REP_ISOFORM')) %>% 
            dplyr::select(!c(n_cells_to_sample)) %>%
            dplyr::filter(Intensity_MeanIntensity_MCHERRY_BACKCORR <= CHERRY_threshold) %>% 
            ungroup() %>% 
            group_by(SET_PAIR_REP_ISOFORM) %>% 
            nest() %>%            
            ungroup() %>% 
            mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% NEG3$SET_PAIR_REP_ISOFORM, NEG3$n_cells_to_sample, NA)) %>%
            mutate(samp = map2(data, n, sample_n)) %>% 
            select(-data) %>%
            unnest(samp) %>% 
            dplyr::select(!n),silent = TRUE)
        
        if(inherits(neg_a, "try-error") | iso_cherry_match == FALSE){
          NEG4 <- NEG_DF %>% 
            ungroup() %>% 
            group_by(SET_PAIR_REP_ISOFORM) %>% 
            nest() %>%            
            ungroup() %>% 
            mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% NEG3$SET_PAIR_REP_ISOFORM, NEG3$n_cells_to_sample, NA)) %>%
            mutate(samp = map2(data, n, sample_n)) %>% 
            select(-data) %>%
            unnest(samp) %>% 
            dplyr::select(!n)
          neg_cherry_match <- FALSE
        }
        if(inherits(neg_a, 'tbl') & iso_cherry_match == TRUE){
          NEG4 <- neg_a
          neg_cherry_match <- TRUE
        }
      }
      if(MCHERRY_thresholding_with_max_mean_when_downsampling == FALSE){
        NEG4 <- NEG_DF %>% 
          ungroup() %>% 
          group_by(SET_PAIR_REP_ISOFORM) %>% 
          nest() %>%            
          ungroup() %>% 
          mutate(n = if_else(SET_PAIR_REP_ISOFORM %in% NEG3$SET_PAIR_REP_ISOFORM, NEG3$n_cells_to_sample, NA)) %>%
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data) %>%
          unnest(samp) %>% 
          dplyr::select(!n)
        neg_cherry_match <- FALSE
      }
      #Negative control Count
      n_negative <- NEG4 %>% 
        ungroup() %>% 
        group_by(SET_PAIR) %>% 
        count() %>% 
        dplyr::pull(n)
      
      #
      NEGATIVE_DATASET_TO_TRIM <- proportional_sampling(X =  NEG4,
                                                        Y = n_postrans_for_balanced_control, # Relies on the isoform one
                                                        seed = NULL)
      
    }
  }
    #Positive control
    if(include_positive_control == TRUE){
      POSITIVE_DATASET_TO_TRIM <- proportional_sampling(X =  POS4,
                                                        Y = n_postrans_for_balanced_control, # Relies on the isoform one
                                                        seed = NULL)
    }
  
  ###
  if(downsample_for_equal_classes == FALSE && sample_controls_to_match_pair == TRUE){
    n_iso <- DATASET %>% 
      dplyr::mutate(GENE_PAIR = word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))) %>% 
      dplyr::filter(GENE_PAIR == PAIR) %>%
      dim()
    n_iso <- n_iso[1]
    
    n_negative <- DATASET %>% 
      dplyr::mutate(GENE_PAIR = word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))) %>% 
      dplyr::filter(GENE_PAIR == 'PAIR_RAB5A') %>%
      dim()
    n_negative <- n_negative[1]
    if(n_negative > n_iso){
      NEGATIVE_DATASET_TO_TRIM <- DATASET %>%
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A') %>% 
        slice_sample(n = n_iso)
    }
    if(n_negative <= n_iso){
      NEGATIVE_DATASET_TO_TRIM <- DATASET %>% 
        ungroup() %>%
        dplyr::filter(word(SET_PAIR, start = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A')
    }
    if(include_positive_control == TRUE){
      n_positive <- DATASET %>%
        ungroup() %>% 
        dplyr::filter(word(SET_PAIR, star = 3, end = 4, sep = fixed('_')) == 'PAIR_RAB5A' | word(SET_PAIR, star = 3, end = 4, sep = fixed('_')) == 'PAIR_MCHERRY') %>% 
        dplyr::mutate(ISOFORM = case_when(
          ISOFORM %in% c('WT','DEx') ~ 'RAB5A',
          ISOFORM == 'MCHERRY' ~ 'MCHERRY'),
          SET_PAIR = paste0(Experiment_Set,'_PAIR_POSITIVE'),
          SET_PAIR_REP_ISOFORM = paste0(SET_PAIR,'_',Experiment_Rep,'_',ISOFORM)) %>% 
        dim()
      n_positive <- n_positive[1]
      
      if(n_positive > n_iso){
        POSITIVE_DATASET_TO_TRIM <- POSITIVE_DF %>% 
          slice_sample(n = n_iso)
      }
      if(n_positive <= n_iso){
        POSITIVE_DATASET_TO_TRIM <- POSITIVE_DF
      }
    }
  }
  #####
  ### Downsampling positive control (currently using RAB5A cells from both wells)
  #####
  
  #####
  ### Creating variable with correct name if downsampling was not done
  if(downsample_for_equal_classes == FALSE){
    ISO_DATASET_TO_TRIM <- DATASET
  }
  ###
  ISOFORM_DATASET <- ISO_DATASET_TO_TRIM %>% 
    dplyr::mutate(GENE_PAIR = word(SET_PAIR, start = 3, end = 4, sep = fixed("_"))) %>% 
    dplyr::filter(GENE_PAIR == PAIR) %>%
    dplyr::select(ISOFORM, all_of(FEATURE_START:length(colnames(DATASET)))) %>%
    dplyr::select(!c(starts_with('Intensity') & ends_with('MCHERRY_BACKCORR') & !contains('MeanIntensity_'))) #Removes intensity data (although leaves intensities normalised to mean intensities)
  #
  if(include_negative_control == TRUE){
    NEGATIVE_DATASET <- NEGATIVE_DATASET_TO_TRIM %>% 
      dplyr::select(ISOFORM, all_of(FEATURE_START:length(colnames(DATASET)))) %>%
      dplyr::select(!c(starts_with('Intensity') & ends_with('MCHERRY_BACKCORR') & !contains('MeanIntensity_'))) #Removes intensity data (although leaves intensities normalised to mean intensities)
  }
  # 
  if(include_positive_control == TRUE){
    POSITIVE_DATASET <- POSITIVE_DATASET_TO_TRIM %>% 
      dplyr::select(ISOFORM, all_of(FEATURE_START:length(colnames(DATASET)))) %>%
      dplyr::select(!c(starts_with('Intensity') & ends_with('MCHERRY_BACKCORR') & !contains('MeanIntensity_'))) #Removes intensity data (although leaves intensities normalised to mean intensities)
  }
  #
  if(include_positive_control == FALSE & include_negative_control == FALSE){
    dataset_list <- list('ISOFORM_DF' = ISOFORM_DATASET)
  }
  #
  if(include_positive_control == FALSE & include_negative_control == TRUE){
    dataset_list <- list('ISOFORM_DF' = ISOFORM_DATASET, 'NEGATIVE_DF' = NEGATIVE_DATASET)
  }
  #
  if(include_positive_control == TRUE & include_negative_control == TRUE){
    dataset_list <- list('ISOFORM_DF' = ISOFORM_DATASET, 'NEGATIVE_DF' = NEGATIVE_DATASET, 'POSITIVE_DF' = POSITIVE_DATASET)
  }
  #return(dataset_list)

  #Caret models
  #####
  #Control object for caret models.   
  if(length(unique(ISOFORM_DATASET$ISOFORM)) == 2){
    myControl <- trainControl(
      method = "cv", #cv = cross validation
      number = 10, #Number of folds
      summaryFunction = twoClassSummary,
      classProbs = TRUE, # <- Super important!
      verboseIter = print_training_progress,
      preProcOptions = list(thresh = 0.9, freqCut = 70/30, uniqueCut = 25), #Note, this is the way to modify how the arguments in the proprocess function works when used within the train workflow
    )
  }
  if(length(unique(ISOFORM_DATASET$ISOFORM)) > 2){
    myControl <- trainControl(
      method = "cv", #cv = cross validation
      number = 10, #Number of folds
      summaryFunction = defaultSummary,
      classProbs = TRUE, # <- Super important!
      verboseIter = print_training_progress,
      preProcOptions = list(thresh = 0.9, freqCut = 70/30, uniqueCut = 25), #Note, this is the way to modify how the arguments in the proprocess function works when used within the train workflow
    )
  }
  
  if(model_type == 'glmnet'){
    #Training glmnet model
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      preProcess = c('zv','nzv','center','scale'),
      trControl = myControl
    )
    #  negative_model <- train(
    #    ISOFORM ~ .,
    #    data = NEGATIVE_DATASET,
    #    method = model_type,
    #    preProcess = c('zv','nzv','center','scale'),
    #    trControl = myControl
    #  )
    #  if(include_positive_control == TRUE){
    #    positive_model <- train(
    #      ISOFORM ~ .,
    #      data = POSITIVE_DATASET,
    #      method = model_type,
    #      preProcess = c('zv','nzv','center','scale'),
    #      trControl = myControl
    #    )
    # }
  }
  #Svm classifiers - can be extended if you want other kernels etc.
  if(model_type %in% c('svmLinear2')){
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      preProcess = c('zv','center','scale','pca'),
      trControl = myControl
    )
  }
  #Ranger (a random-forest type)
  if(model_type %in% c('ranger','rf')){
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      preProcess = c('zv','center','scale'),
      trControl = myControl
    )
    negative_model <- train(
      ISOFORM ~ .,
      data = NEGATIVE_DATASET,
      method = model_type,
      preProcess = c('zv','center','scale'),
      trControl = myControl
    )
    if(include_positive_control == TRUE){
      positive_model <- train(
        ISOFORM ~ .,
        data = POSITIVE_DATASET,
        method = model_type,
        preProcess = c('zv','center','scale'),
        trControl = myControl
      )
    }
  }
  #Neural network (feedforward neural network, FNN)
  if(model_type == 'nnet'){
    isoform_model <- train(
      ISOFORM ~ .,
      data = ISOFORM_DATASET,
      method = model_type,
      MaxNWts = 1060,
      preProcess = c('zv','center','scale'),
      trControl = myControl
    )
    if(include_negative_control == TRUE){
      negative_model <- train(
        ISOFORM ~ .,
        data = NEGATIVE_DATASET,
        method = model_type,
        MaxNWts = 1060,
        preProcess = c('zv','center','scale'),
        trControl = myControl
      )
    }
    if(include_positive_control == TRUE){
      positive_model <- train(
        ISOFORM ~ .,
        data = POSITIVE_DATASET,
        method = model_type,
        MaxNWts = 1060,
        preProcess = c('zv','center','scale'),
        trControl = myControl
      )
    }
  }
  #####
  if(reduce_model_output == TRUE){
    reduced_model_lists <- function(LIST){
      #
      variable_importance <- varImp(LIST)$importance
      #
      LIST$trainingData <- LIST$trainingData %>% 
        dplyr::select(.outcome, Intensity_MeanIntensity_MCHERRY_BACKCORR)
      #
      reduced_list <- list('results' = LIST$results, 'training_data_intensity' = LIST$trainingData, 'variable_importance' = variable_importance)
      return(reduced_list)
    }
    #
    isoform_model <- reduced_model_lists(isoform_model)
    isoform_model <- append(isoform_model, values = c(list('iso_cherry_match' = iso_cherry_match),list('SET_PAIR' = iso_set_pair)))
    
    #
    if(include_negative_control == TRUE){
      negative_model <- reduced_model_lists(negative_model)
      negative_model <- append(negative_model, values = c(list('neg_cherry_match' = neg_cherry_match),list('SET_PAIR' = neg_set_pair)))
    }
    #
    if(include_positive_control == TRUE){
      positive_model <- reduced_model_lists(positive_model)
      positive_model <- append(positive_model, values = c(list('pos_cherry_match' = pos_cherry_match),list('SET_PAIR' = pos_set_pair)))
    }
  }
  if(include_positive_control == FALSE & include_negative_control == FALSE){
    model_list <- list('ISOFORM_MODEL' = isoform_model)
  }
  if(include_positive_control == FALSE & include_negative_control == TRUE){
    model_list <- list('ISOFORM_MODEL' = isoform_model, 'NEGATIVE_MODEL' = negative_model)
  }
  if(include_positive_control == TRUE & include_negative_control == TRUE){
    model_list <- list('ISOFORM_MODEL' = isoform_model, 'NEGATIVE_MODEL' = negative_model, 'POSITIVE_MODEL' = positive_model)
  }
  return(model_list)
}

# Replications of sampling + classifier training + assessment
#Pairs to be test for the various sets
#Smaller sets of pairs for set 1 and A, separately. These are used for the positive classifiers
PAIRS_1 <- c('PAIR_RAB5A', 'PAIR_PosControl','PAIR_01','PAIR_03','PAIR_04','PAIR_08','PAIR_09')
PAIRS_A <- c('PAIR_RAB5A', 'PAIR_PosControl','PAIR_01','PAIR_03','PAIR_04','PAIR_08','PAIR_09')
#These are the actual isoform tests
PAIRS_1A <- c('PAIR_RAB5A', 'PAIR_PosControl','PAIR_01','PAIR_02','PAIR_03','PAIR_04','PAIR_05','PAIR_06','PAIR_07','PAIR_08','PAIR_09','PAIR_10')
PAIRS_1A_alternative_comps <- c('PAIR_03DExVSPanEx','PAIR_03bothEx','PAIR_03ExVsEx','PAIR_03Ex1VsBoth')
PAIRS_B <- c('PAIR_01','PAIR_02','PAIR_03','PAIR_04','PAIR_05','PAIR_06','PAIR_07','PAIR_08','PAIR_09','PAIR_10','PAIR_RAB5A','PAIR_PosControl')
PAIRS_C <- c('PAIR_02','PAIR_04','PAIR_05','PAIR_06','PAIR_07','PAIR_08','PAIR_09','PAIR_10','PAIR_RAB5A', 'PAIR_PosControl')
#
PAIRS_D_wo_01 <- c('PAIR_RAB5A', 'PAIR_PosControl','PAIR_02','PAIR_03','PAIR_04','PAIR_05','PAIR_06','PAIR_07','PAIR_08','PAIR_09')#
PAIRS_D_alternative_comps <- c('PAIR_01wtVS57','PAIR_01wtVS168','PAIR_01dex57VS168','PAIR_08DExVSaltEx','PAIR_08bothEx','PAIR_08ExVsEx','PAIR_08Ex1VsBoth')


#
repeat_sampling_and_classifier <- function(CELL_DATASET, PAIRS_LIST, SET_ID, NUMBER_OF_REPLICATIONS, SAVE_DIRECTORY){
  for(i in PAIRS_LIST){
    set.seed(123) #To make it consistent upon re-tries
    replicate_nnet <- replicate(n = NUMBER_OF_REPLICATIONS, expr = isoform_classify(DATASET = CELL_DATASET,
                                                                                    PAIR = i,
                                                                                    FEATURE_START = 14,
                                                                                    model_type = 'nnet',
                                                                                    downsample_for_equal_classes = TRUE,
                                                                                    sample_controls_to_match_pair = TRUE,
                                                                                    include_positive_control = FALSE,
                                                                                    MCHERRY_thresholding_with_max_mean_when_downsampling = TRUE,
                                                                                    reduce_model_output = TRUE,
                                                                                    print_training_progress = FALSE,
                                                                                    include_negative_control = FALSE))
    assign(paste0(SET_ID,'_',i,'_model_results'),replicate_nnet)
    print(paste0('Saving models for ',i))
    save(replicate_nnet, file = paste0(SAVE_DIRECTORY,'/',SET_ID,'_',i))
    time_now <- Sys.time()
    print(paste0('Starting next pair. Time is ',time_now))
    gc()
  }
}

#Repeat sampling + classifier for NegativeTransfection dataset
repeat_sampling_and_classifier_NegTransfected <- function(NEGTRANS_DATASET, POSTRANS_DATASET, PAIRS_LIST, SET_ID, NUMBER_OF_REPLICATIONS, SAVE_DIRECTORY){
  for(i in PAIRS_LIST){
    set.seed(123) #To make it consistent upon re-tries
    replicate_nnet <- replicate(n = NUMBER_OF_REPLICATIONS, expr = isoform_classify_for_negative_cells(DATASET = NEGTRANS_DATASET,
                                                                                    PAIR = i,
                                                                                    FEATURE_START = 14,
                                                                                    model_type = 'nnet',
                                                                                    downsample_for_equal_classes = TRUE,
                                                                                    sample_controls_to_match_pair = TRUE,
                                                                                    include_positive_control = FALSE,
                                                                                    MCHERRY_thresholding_with_max_mean_when_downsampling = TRUE,
                                                                                    reduce_model_output = TRUE,
                                                                                    print_training_progress = FALSE,
                                                                                    POSITIVE_DATASET = POSTRANS_DATASET,
                                                                                    include_negative_control = FALSE))
    assign(paste0('NegTransfected_',SET_ID,'_',i,'_model_results'),replicate_nnet)
    print(paste0('Saving models for ',i))
    save(replicate_nnet, file = paste0(SAVE_DIRECTORY,'/NegTransModel_',SET_ID,'_',i))
    time_now <- Sys.time()
    print(paste0('Starting next pair. Time is ',time_now))
    gc()
  }
}

#Positive transfectants models
#####
#Define directory to save models
Directory_to_save_postransfected_models_approach1_lowerbin <- 'INSERT DIRECTORY'
Directory_to_save_postransfected_models_approach1_upperbin <- 'INSERT DIRECTORY'
#
Directory_to_save_postransfected_models_approach2_lowerbin <- 'INSERT DIRECTORY'
Directory_to_save_postransfected_models_approach2_upperbin <- 'INSERT DIRECTORY'
#
Directory_to_save_postransfected_models_approach3_lowerbin <- 'INSERT DIRECTORY'
Directory_to_save_postransfected_models_approach3_upperbin <- 'INSERT DIRECTORY'

#####

#SET 1A, approach 1
#Lower
#repeat_sampling_and_classifier(CELL_DATASET = ,
#                               PAIRS_LIST = PAIRS_1,
#                               SET_ID = 'SET_1',
#                               NUMBER_OF_REPLICATIONS = 50,
#                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
#repeat_sampling_and_classifier(CELL_DATASET = ,
#                               PAIRS_LIST = PAIRS_1,
#                               SET_ID = 'SET_1',
#                               NUMBER_OF_REPLICATIONS = 50,
#                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)


### SET 1A, approach 1 ###
#####
#SET 1A, approach 1 
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_ap1_lowerbin,
                               PAIRS_LIST = PAIRS_1A,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_ap1_upperbin,
                               PAIRS_LIST = PAIRS_1A,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)
#SET 1A, approach 2
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_ap2_lowerbin,
                               PAIRS_LIST = PAIRS_1A,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_ap2_upperbin,
                               PAIRS_LIST = PAIRS_1A,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_upperbin)
#SET 1A, approach 3
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_ap3_lowerbin,
                               PAIRS_LIST = PAIRS_1A,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_ap3_upperbin,
                               PAIRS_LIST = PAIRS_1A,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_upperbin)

#####

### Set 1A, alternative comparisons ###
#####
## SET 1A alt comps, approach 1 ##
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_altcomps_ap1_lowerbin,
                               PAIRS_LIST = PAIRS_1A_alternative_comps,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_altcomps_ap1_upperbin,
                               PAIRS_LIST = PAIRS_1A_alternative_comps,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)
## SET 1A alt comps, approach 2 ##
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_altcomps_ap2_lowerbin,
                               PAIRS_LIST = PAIRS_1A_alternative_comps,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_altcomps_ap2_upperbin,
                               PAIRS_LIST = PAIRS_1A_alternative_comps,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_upperbin)
## SET 1A alt comps, approach 3 ##
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_altcomps_ap3_lowerbin,
                               PAIRS_LIST = PAIRS_1A_alternative_comps,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_1A_altcomps_ap3_upperbin,
                               PAIRS_LIST = PAIRS_1A_alternative_comps,
                               SET_ID = 'SET_1A',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_upperbin)
#####

### Set B ###
#####
#SET B, approach 1 
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_B_ap1_lowerbin,
                               PAIRS_LIST = PAIRS_B,
                               SET_ID = 'SET_B',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_B_ap1_upperbin,
                               PAIRS_LIST = PAIRS_B,
                               SET_ID = 'SET_B',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)
#SET B, approach 2 
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_B_ap2_lowerbin,
                               PAIRS_LIST = PAIRS_B,
                               SET_ID = 'SET_B',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_B_ap2_upperbin,
                               PAIRS_LIST = PAIRS_B,
                               SET_ID = 'SET_B',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_upperbin)
#SET B, approach 3 
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_B_ap3_lowerbin,
                               PAIRS_LIST = PAIRS_B,
                               SET_ID = 'SET_B',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_B_ap3_upperbin,
                               PAIRS_LIST = PAIRS_B,
                               SET_ID = 'SET_B',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_upperbin)
#####

### Set C ###
#####
#SET C, approach 1
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_C_ap1_lowerbin,
                               PAIRS_LIST = PAIRS_C,
                               SET_ID = 'SET_C',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_C_ap1_upperbin,
                               PAIRS_LIST = PAIRS_C,
                               SET_ID = 'SET_C',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)

#SET C, approach 2
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_C_ap2_lowerbin,
                               PAIRS_LIST = PAIRS_C,
                               SET_ID = 'SET_C',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_C_ap2_upperbin,
                               PAIRS_LIST = PAIRS_C,
                               SET_ID = 'SET_C',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_upperbin)
#SET C, approach 3
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_C_ap3_lowerbin,
                               PAIRS_LIST = PAIRS_C,
                               SET_ID = 'SET_C',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_C_ap3_upperbin,
                               PAIRS_LIST = PAIRS_C,
                               SET_ID = 'SET_C',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_upperbin)
#####

### SET D ###
#####
# Set D, approach 1
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_ap1_lowerbin,
                               PAIRS_LIST = PAIRS_D_wo_01,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_ap1_upperbin,
                               PAIRS_LIST = PAIRS_D_wo_01,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)
# Set D, approach 2
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_ap2_lowerbin,
                               PAIRS_LIST = PAIRS_D_wo_01,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_ap2_upperbin,
                               PAIRS_LIST = PAIRS_D_wo_01,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_upperbin)
# Set D, approach 3
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_ap3_lowerbin,
                               PAIRS_LIST = PAIRS_D_wo_01,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_ap3_upperbin,
                               PAIRS_LIST = PAIRS_D_wo_01,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_upperbin)
#####

### SET D, ALTERNATIVE COMPARISONS ###
#####
# Set D, approach 1
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_altcomps_ap1_lowerbin,
                               PAIRS_LIST = PAIRS_D_alternative_comps,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_altcomps_ap1_upperbin,
                               PAIRS_LIST = PAIRS_D_alternative_comps,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach1_upperbin)
# Set D, approach 2
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_altcomps_ap2_lowerbin,
                               PAIRS_LIST = PAIRS_D_alternative_comps,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_altcomps_ap2_upperbin,
                               PAIRS_LIST = PAIRS_D_alternative_comps,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach2_upperbin)
# Set D, approach 3
#Lower
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_altcomps_ap3_lowerbin,
                               PAIRS_LIST = PAIRS_D_alternative_comps,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier(CELL_DATASET = Transfected_for_class_D_altcomps_ap3_upperbin,
                               PAIRS_LIST = PAIRS_D_alternative_comps,
                               SET_ID = 'SET_D',
                               NUMBER_OF_REPLICATIONS = 50,
                               SAVE_DIRECTORY = Directory_to_save_postransfected_models_approach3_upperbin)
#####


#Negative transfectants models
#####
#Define directory to save models
Directory_to_save_negtransfected_models_approach1_lowerbin <- 'INSERT DIRECTORY'
Directory_to_save_negtransfected_models_approach1_upperbin <- 'INSERT DIRECTORY'
#
Directory_to_save_negtransfected_models_approach2_lowerbin <- 'INSERT DIRECTORY'
Directory_to_save_negtransfected_models_approach2_upperbin <- 'INSERT DIRECTORY'
#
Directory_to_save_negtransfected_models_approach3_lowerbin <- 'INSERT DIRECTORY'
Directory_to_save_negtransfected_models_approach3_upperbin <- 'INSERT DIRECTORY'

### SET 1A ###
#####
# SET 1A, approach 1
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A,
                                              POSTRANS_DATASET = Transfected_for_class_1A_ap1_lowerbin,
                                              PAIRS_LIST = PAIRS_1A,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A,
                                              POSTRANS_DATASET = Transfected_for_class_1A_ap1_upperbin,
                                              PAIRS_LIST = PAIRS_1A,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_upperbin)
#SET 1A, approach 2
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A,
                                              POSTRANS_DATASET = Transfected_for_class_1A_ap2_lowerbin,
                                              PAIRS_LIST = PAIRS_1A,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A,
                                              POSTRANS_DATASET = Transfected_for_class_1A_ap2_upperbin,
                                              PAIRS_LIST = PAIRS_1A,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_upperbin)
#SET 1A, approach 3
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A,
                                              POSTRANS_DATASET = Transfected_for_class_1A_ap3_lowerbin,
                                              PAIRS_LIST = PAIRS_1A,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A,
                                              POSTRANS_DATASET = Transfected_for_class_1A_ap3_upperbin,
                                              PAIRS_LIST = PAIRS_1A,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_upperbin)
#####

### SET 1A alternative comps ###
#####
#SET 1A alternative comps, approach 1
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_1A_altcomps_ap1_lowerbin,
                                              PAIRS_LIST = PAIRS_1A_alternative_comps,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_1A_altcomps_ap1_upperbin,
                                              PAIRS_LIST = PAIRS_1A_alternative_comps,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_upperbin)
#SET 1A alternative comps, approach 2
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_1A_altcomps_ap2_lowerbin,
                                              PAIRS_LIST = PAIRS_1A_alternative_comps,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_1A_altcomps_ap2_upperbin,
                                              PAIRS_LIST = PAIRS_1A_alternative_comps,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_upperbin)
#SET 1A alternative comps, approach 3
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_1A_altcomps_ap3_lowerbin,
                                              PAIRS_LIST = PAIRS_1A_alternative_comps,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_1A_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_1A_altcomps_ap3_upperbin,
                                              PAIRS_LIST = PAIRS_1A_alternative_comps,
                                              SET_ID = 'SET_1A',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_upperbin)
#####




### SET B ###
#####
#SET B, approach 1
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_B,
                                              POSTRANS_DATASET = Transfected_for_class_B_ap1_lowerbin,
                                              PAIRS_LIST = PAIRS_B,
                                              SET_ID = 'SET_B',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_B,
                                              POSTRANS_DATASET = Transfected_for_class_B_ap1_upperbin,
                                              PAIRS_LIST = PAIRS_B,
                                              SET_ID = 'SET_B',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_upperbin)
#SET B, approach 2
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_B,
                                              POSTRANS_DATASET = Transfected_for_class_B_ap2_lowerbin,
                                              PAIRS_LIST = PAIRS_B,
                                              SET_ID = 'SET_B',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_B,
                                              POSTRANS_DATASET = Transfected_for_class_B_ap2_upperbin,
                                              PAIRS_LIST = PAIRS_B,
                                              SET_ID = 'SET_B',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_upperbin)
#SET B, approach 3
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_B,
                                              POSTRANS_DATASET = Transfected_for_class_B_ap3_lowerbin,
                                              PAIRS_LIST = PAIRS_B,
                                              SET_ID = 'SET_B',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_lowerbin)

#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_B,
                                              POSTRANS_DATASET = Transfected_for_class_B_ap3_upperbin,
                                              PAIRS_LIST = PAIRS_B,
                                              SET_ID = 'SET_B',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_upperbin)
#####


### SET C ###
#####
#SET C, approach 1
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_C,
                                              POSTRANS_DATASET = Transfected_for_class_C_ap1_lowerbin,
                                              PAIRS_LIST = PAIRS_C,
                                              SET_ID = 'SET_C',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_C,
                                              POSTRANS_DATASET = Transfected_for_class_C_ap1_upperbin,
                                              PAIRS_LIST = PAIRS_C,
                                              SET_ID = 'SET_C',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_upperbin)
#SET C, approach 2
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_C,
                                              POSTRANS_DATASET = Transfected_for_class_C_ap2_lowerbin,
                                              PAIRS_LIST = PAIRS_C,
                                              SET_ID = 'SET_C',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_C,
                                              POSTRANS_DATASET = Transfected_for_class_C_ap2_upperbin,
                                              PAIRS_LIST = PAIRS_C,
                                              SET_ID = 'SET_C',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_upperbin)
#SET C, approach 3
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_C,
                                              POSTRANS_DATASET = Transfected_for_class_C_ap3_lowerbin,
                                              PAIRS_LIST = PAIRS_C,
                                              SET_ID = 'SET_C',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_C,
                                              POSTRANS_DATASET = Transfected_for_class_C_ap3_upperbin,
                                              PAIRS_LIST = PAIRS_C,
                                              SET_ID = 'SET_C',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_upperbin)
#####


### SET D ###
#####
#SET D, approach 1
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D,
                                              POSTRANS_DATASET = Transfected_for_class_D_ap1_lowerbin,
                                              PAIRS_LIST = PAIRS_D_wo_01,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D,
                                              POSTRANS_DATASET = Transfected_for_class_D_ap1_upperbin,
                                              PAIRS_LIST = PAIRS_D_wo_01,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_upperbin)
#SET D, approach 2
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D,
                                              POSTRANS_DATASET = Transfected_for_class_D_ap2_lowerbin,
                                              PAIRS_LIST = PAIRS_D_wo_01,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D,
                                              POSTRANS_DATASET = Transfected_for_class_D_ap2_upperbin,
                                              PAIRS_LIST = PAIRS_D_wo_01,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_upperbin)
#SET D, approach 3
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D,
                                              POSTRANS_DATASET = Transfected_for_class_D_ap3_lowerbin,
                                              PAIRS_LIST = PAIRS_D_wo_01,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D,
                                              POSTRANS_DATASET = Transfected_for_class_D_ap3_upperbin,
                                              PAIRS_LIST = PAIRS_D_wo_01,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_upperbin)

#####

### SET D ALTERNATIVE COMPARISONS ###
#SET D ALTERNATIVE COMPS, approach 1
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_D_altcomps_ap1_lowerbin,
                                              PAIRS_LIST = PAIRS_D_alternative_comps,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_D_altcomps_ap1_upperbin,
                                              PAIRS_LIST = PAIRS_D_alternative_comps,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach1_upperbin)
#SET D ALTERNATIVE COMPS, approach 2
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_D_altcomps_ap2_lowerbin,
                                              PAIRS_LIST = PAIRS_D_alternative_comps,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_D_altcomps_ap2_upperbin,
                                              PAIRS_LIST = PAIRS_D_alternative_comps,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach2_upperbin)
#SET D ALTERNATIVE COMPS, approach 3
#Lower
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_D_altcomps_ap3_lowerbin,
                                              PAIRS_LIST = PAIRS_D_alternative_comps,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_lowerbin)
#Upper
repeat_sampling_and_classifier_NegTransfected(NEGTRANS_DATASET = Negtrans_for_class_D_altcomp,
                                              POSTRANS_DATASET = Transfected_for_class_D_altcomps_ap3_upperbin,
                                              PAIRS_LIST = PAIRS_D_alternative_comps,
                                              SET_ID = 'SET_D',
                                              NUMBER_OF_REPLICATIONS = 50,
                                              SAVE_DIRECTORY = Directory_to_save_negtransfected_models_approach3_upperbin)


