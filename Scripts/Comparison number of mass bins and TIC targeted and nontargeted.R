#############################################################################################################################################
library(tidyverse)
#load required datafile
load(file = "./RDA/PseudoMALDI_targeted_extended.rda")
load(file= "./RDA/SNOG_complete_diagnostic.rda")
load(file= "./RDA/SNOG_Precursor_gather.rda")
load(file= "./RDA/Joined_TIC.rda")

##Comparison of mass bins targeted vs nontargeted, tissue dependent


#Number of unique masses per tissue in targeted data
n_targeted_tissue_dependent <- gather(PseudoMALDI_targeted_extended, key = "Sample", value = "intensity", - "mass") %>% mutate(mass = round(mass, digits = 1)) %>%
  filter(intensity > 10^7, mass < 4000, mass > 1000) %>% group_by(Sample, mass) %>% slice(1) %>%  #insert intensity threshold here
  ungroup() %>% group_by(Sample) %>% summarize(n_targeted = n()) %>%
  mutate(Sample = tolower(Sample)) %>% mutate(Sample = recode(Sample, "mammaryggland2" = "mammarygland2"))

#Filter nontargeted data for SNOG score > 5; do not remove mutiply occuring masses yet
SNOG_TIC_n <- as.data.frame(SNOG_complete_diagnostic %>%
                              mutate(SNOG = log10(((X224.11177 -             ##calculate SNOGscore
                                                      ((X183.08629+X167.09139+X312.12889+X328.12379)*10))/BPI)*TIC)) %>% 
                              select(-SNOGscore) %>% 
                              mutate(Sample = recode(Sample, "BrownAdiposeTissue1" = "BAT1", "BrownAdiposeTissue2" = "BAT2", #Recode names for joining
                                                     "colon1re" = "colon1", "WhiteAdiposeTissue1" = "WAT1", "WhiteAdiposeTissue2" = "WAT2")) %>%
                              filter(SNOG > 5) %>% 
                              group_by(Sample, M))

#get precursor intensity for the filtered MSMS spectra per tissue
j <- inner_join(SNOG_TIC_n,SNOG_Precursor_gather, by = c("M", "Sample") )

#filter for desired intensity value BEFORE removing mutiply occuring masses, count the number of unique mass bins per tissue
j_filter <- j %>% filter(intensity > 10^7) %>% group_by(Sample, M) %>% do(head(.,1)) %>% ungroup() %>%    #insert intensity threshold here
  group_by(Sample) %>% filter(M < 4000, M > 1000) %>% summarize(n_nontargeted = n()) %>%
  mutate(Sample = tolower(Sample))

#join the mass bins from targeted and untargeted approach into 1 data frame
comparison <- inner_join(n_targeted_tissue_dependent, j_filter, by = "Sample")

#Plot the data
comparison %>% ggplot() + geom_col(aes(Sample, n_targeted), fill = "blue", alpha = 0.7) + 
  geom_col(aes(Sample, n_nontargeted), fill = "orange", alpha = 0.7) + theme_classic() +
  xlab("") + ylab("number of unique mass-bins") + 
  theme(text = element_text(size=16), axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5)) 


######################################################################################################################################
###########comparison TIC targetetd vs non-targeted tissue dependent


#Summed intensities of targeted data, stratified by tissue
targeted_summed <- PseudoMALDI_targeted_extended %>% gather(key = "Sample", value = "intensity", - mass) %>% filter(intensity > 0) %>%   #insert intensity threshold here
  group_by(Sample) %>% summarize(TIC_targeted = sum(intensity))  %>% mutate(Sample = tolower(Sample)) %>%
  mutate(Sample = recode(Sample, "mammaryggland2" = "mammarygland2")) 


#nontargeted data
SNOG_TIC <- as.data.frame(SNOG_complete_diagnostic %>%
                            mutate(SNOG = log10(((X224.11177 -             
                                                    ((X183.08629+X167.09139+X312.12889+X328.12379)*10))/BPI)*TIC)) %>% 
                            select(-SNOGscore) %>% 
                            mutate(Sample = recode(Sample, "BrownAdiposeTissue1" = "BAT1", "BrownAdiposeTissue2" = "BAT2",
                                                   "colon1re" = "colon1", "WhiteAdiposeTissue1" = "WAT1", "WhiteAdiposeTissue2" = "WAT2")) %>%
                            filter(SNOG > 5) %>% 
                            group_by(Sample, M) %>%
                            do(head(.,1))) %>% filter(M < 4000, M > 1000)

Joined_TIC <- inner_join(SNOG_Precursor_gather, SNOG_TIC, by = c("Sample", "M")) 

#Sum of the intensities of the nontargeted data, combined into data frame with summed intensties from targeted data
Joined_TIC_sum <- Joined_TIC %>% filter(intensity > 0) %>% group_by(Sample) %>% summarize(TIC = sum(intensity)) %>% #insert intensity threshold here
  mutate(Sample = tolower(Sample)) %>% inner_join(., targeted_summed, by = "Sample")

#Visualization of summed intensities from nontargeted and targeted data
Joined_TIC_correlation_sum %>% ggplot()  + geom_col(aes(Sample,TIC_targeted), fill = "blue", alpha = 0.7) +
  geom_col(aes(Sample,TIC), fill = "orange", alpha = 0.7)+   theme_classic() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5)) + xlab("") + ylab("TIC")



#######################################################################################################################
#####Getting the numbers for Venn diagrams

#Number of unique masses  in targeted data
targeted <- gather(PseudoMALDI_targeted_extended, key = "Sample", value = "intensity", - "mass") %>% mutate(mass = round(mass, digits = 1)) %>%
  filter(intensity > 0, mass < 4000, mass > 1000) %>% group_by(mass) %>% slice(1) %>%  #insert intensity threshold here
  ungroup() %>% mutate(M = mass)

targeted_threshold <- gather(PseudoMALDI_targeted_extended, key = "Sample", value = "intensity", - "mass") %>% mutate(mass = round(mass, digits = 1)) %>%
  filter(intensity > 10^7, mass < 4000, mass > 1000) %>% group_by(mass) %>% slice(1) %>%  #insert intensity threshold here
  ungroup() %>% mutate(M = mass)


#TIC untargeted
sum(Joined_TIC$intensity) # without intensity filter
sum(Joined_TIC_threshold$intensity) #with intensity filter

#TIC matched
targeted_list <- targeted %>% select(M)
TIC_matched <- inner_join(targeted_list, Joined_TIC, by = "M")
sum(TIC_matched$intensity) #TIC matched without intensity filter

targeted_list_threshold <- targeted_threshold %>% select(M)
TIC_matched_threshold <- inner_join(targeted_list_threshold, Joined_TIC_threshold, by = "M")
sum(TIC_matched_threshold$intensity) #TIC matched with intensity filter 10^7

#TIC not matched nontargeted
no_match <- anti_join(Joined_TIC,targeted_list , by = "M")
sum(no_match$intensity) #without intensity filter

no_match_threshold <- anti_join(Joined_TIC_threshold,targeted_list_threshold , by = "M")
sum(no_match_threshold$intensity) # with intensity filter 10^7

#TIC not matched targeted
no_match_targeted <- anti_join(targeted, Joined_TIC, by = "M")
sum(no_match_targeted$intensity) #without intensity filter

no_match_targeted_threshold <- anti_join(targeted_threshold, Joined_TIC_threshold, by = "M")
sum(no_match_targeted_threshold$intensity)  #without intensity filter

no_match_targeted %>% group_by(M) %>% do(head(.,1)) %>% nrow()

#targeted mass bins
nrow(targeted) #no threshold
nrow(targeted_threshold) #threshold intensity 10^7

#nontargeted mass bins
Joined_TIC %>% group_by(M) %>% do(head(.,1)) %>% nrow() #no threshold
Joined_TIC_threshold <- Joined_TIC %>% filter(intensity > 10^7)
Joined_TIC_threshold %>% group_by(M) %>% do(head(.,1)) %>% nrow() #threshold

#matched mass bins
matched <- inner_join(targeted, Joined_TIC, by = "M") %>% group_by(M) %>% do(head(.,1)) %>% nrow()
matched_threshold <- inner_join(targeted_threshold, Joined_TIC_threshold, by = "M") %>% group_by(M) %>% do(head(.,1)) %>% nrow()

