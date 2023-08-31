###Data preparation

library(tidyverse)

####Add download URLs####

SNOG_complete_diagnostic <- read.csv("./Files/MS1/SNOG_complete.csv", header = TRUE, sep = ";") #load the data set
SNOG_complete_diagnostic <- SNOG_complete_diagnostic %>% mutate(M = round(M, digits = 1)) #round the precursor mass to make it comparable
PseudoMALDI_targeted_extended <- read.csv("./Files/MS1/PseudoMALDI_targeted_extended.csv") #targeted data set

SNOG_Precursor <- read.csv("./Files/MS1/Oxobased_PseudoMALDI.csv", header = TRUE, sep = ";") # load the data frame with the precursor intensities
SNOG_Precursor_gather <- gather(SNOG_Precursor, key = "Sample", value = intensity, -1) %>%
  mutate(M = round(mass, digits = 1)) %>% select(-mass)  #round to make comparable with data frame to join with

save(SNOG_Precursor_gather, file = "./RDA/SNOG_Precursor_gather.rda") #save as .rda
save(SNOG_complete_diagnostic, file = "./RDA/SNOG_complete_diagnostic.rda") #save as .rda
save(PseudoMALDI_targeted_extended, file = "./RDA/PseudoMALDI_targeted_extended.rda") #save as .rda

####Filter data SNOG-score > 5               
SNOG_TIC <- as.data.frame(SNOG_complete_diagnostic %>%
                            mutate(SNOG = log10(((X224.11177 -      ##calculate SNOGscore
                                                    ((X183.08629+X167.09139+X312.12889+X328.12379)*10))/BPI)*TIC)) %>% 
                            select(-SNOGscore) %>% 
                            mutate(Sample = recode(Sample, "BrownAdiposeTissue1" = "BAT1", "BrownAdiposeTissue2" = "BAT2", #Recode names for joining data frames later
                                                   "colon1re" = "colon1", "WhiteAdiposeTissue1" = "WAT1", "WhiteAdiposeTissue2" = "WAT2")) %>%
                            filter(SNOG > 5) %>% 
                            group_by(M, Sample) %>%
                            do(head(.,1))) %>% filter(M < 4000, M > 1000)  ##if there are multiple MSMS triggered from the same precursor mass, take only the first one per tissue

Joined_TIC <- inner_join(SNOG_Precursor_gather, SNOG_TIC, by = c("Sample", "M")) %>% filter(intensity > 0)#join data frames to get precursor intensities 

save(Joined_TIC, file = "./RDA/Joined_TIC.rda") #save as .rda