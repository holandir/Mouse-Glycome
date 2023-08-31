#load the required packages

library(tidyverse)
library(VennDiagram)

SNOG_complete_diagnostic <- read.csv("./Files/MS1/SNOG_complete.csv", header = TRUE, sep = ";") #load the data set
SNOG_complete_diagnostic <- SNOG_complete_diagnostic %>% mutate(M = round(M, digits = 1)) #round the precursor mass to make it comparable
PseudoMALDI_targeted_extended <- read.csv("./Files/MS1/PseudoMALDI_targeted_extended.csv")

SNOG_Precursor <- read.csv("./Files/MS1/Oxobased_PseudoMALDI.csv", header = TRUE, sep = ";") # load the data frame with the precursor intensities
SNOG_Precursor_gather <- gather(SNOG_Precursor, key = "Sample", value = intensity, -1) %>%
                         mutate(M = round(mass, digits = 1)) %>% select(-mass)  #round to make comparable with data frame to join with

save(SNOG_Precursor_gather, file = "./RDA/SNOG_Precursor_gather.rda") #save as .rda
save(SNOG_complete_diagnostic, file = "./RDA/SNOG_complete_diagnostic.rda") #save as .rda
save(PseudoMALDI_targeted_extended, file = "./RDA/PseudoMALDI_targeted_extended.rda") #save as .rda

####Filter SNOG-score > 5               
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


####################################################################################################
#####Impact of different SNOG-scores on number of unique masses and MSMS spectra# - over all tissues 
load(file ="./RDA/PseudoMALDI_targeted_extended.rda")

#write function for applying different SNOG scores to get unique masses
SNOG_calc <- function(x){
  data <- as.data.frame(SNOG_complete_diagnostic %>%
                          mutate(SNOG = log10(((X224.11177 -            
                          ((X183.08629+X167.09139+X312.12889+X328.12379)*10))/BPI)*TIC)) %>% 
                          select(-SNOGscore) %>% 
                          mutate(Sample = recode(Sample, "BrownAdiposeTissue1" = "BAT1", "BrownAdiposeTissue2" = "BAT2", 
                          "colon1re" = "colon1", "WhiteAdiposeTissue1" = "WAT1", "WhiteAdiposeTissue2" = "WAT2"))) %>%
                          filter(M > 1000, M<4000)
  y <- data %>% filter(x < SNOG) %>% group_by(M) %>% do(head(.,1)) #### x <- SNOG values
  nrow(y) ### count number of mass-bins
}

sequence <- seq(3, 8, 0.2) ##SNOG score values to insert into function
n_SNOG <- sapply(sequence, SNOG_calc) ###apply the different SNOG score values

#write function for applying different SNOG values to get number of MSMS scans
SNOG_calc_n_MSMS <- function(x){
  data <- as.data.frame(SNOG_complete_diagnostic %>%
                          mutate(SNOG = log10(((X224.11177 -            
                          ((X183.08629+X167.09139+X312.12889+X328.12379)*10))/BPI)*TIC)) %>% 
                          select(-SNOGscore) %>% 
                          mutate(Sample = recode(Sample, "BrownAdiposeTissue1" = "BAT1", "BrownAdiposeTissue2" = "BAT2", 
                          "colon1re" = "colon1", "WhiteAdiposeTissue1" = "WAT1", "WhiteAdiposeTissue2" = "WAT2"))) %>%
                          filter(M > 1000, M < 4000)
  y <- data %>% filter(x < SNOG) 
  nrow(y) ### count number of MSMS spectra
}

#apply the different SNOG values (as above)
n_SNOG_MSMS <- sapply(sequence, SNOG_calc_n_MSMS)

###combine the number of MSMS spectra and number of mass bins into one data frame
n_SNOG_complete <- data.frame("SNOG" = sequence, "n" = n_SNOG_MSMS, "n_unique_masses" = n_SNOG )

###plot data
n_SNOG_complete %>% ggplot(aes(x = SNOG)) + geom_point(aes(y = n/26.3), col = "red", alpha = 0.6) + geom_line(aes(y = n/26.3), alpha = 0.2) + 
  geom_point(aes(y = n_unique_masses), col = "blue", alpha = 0.6) + geom_line(aes(y = n_unique_masses), alpha = 0.2) +
  theme_classic() + 
  scale_y_continuous(sec.axis = sec_axis(~ . * 26.3, name = "number of MSMS spectra")) + 
  ylab("number of unique masses") + xlab("SNOG score")
#############################################################################################################################################



#############################################################################################################################################
#Number of masses found with targeted approach
load(file ="./RDA/PseudoMALDI_targeted_extended.rda") #load required datafile

#Number of unique masses per tissue in targeted data
n_targeted_tissue_dependent <- gather(PseudoMALDI_targeted_extended, key = "Sample", value = "intensity", - "mass") %>% mutate(mass = round(mass, digits = 1)) %>%
                               filter(intensity > 10^7, mass < 4000, mass > 1000) %>% group_by(Sample, mass) %>% slice(1) %>%  #insert intensity threshold here
                               ungroup() %>% group_by(Sample) %>% summarize(n_targeted = n()) %>%
                               mutate(Sample = tolower(Sample)) %>% mutate(Sample = recode(Sample, "mammaryggland2" = "mammarygland2"))

#Filter for SNOG score > 5; do not remove mutiply occuring masses yet
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

##########################################################################################################################################


##########################################################################################################################################
#######How many mass bins are there after "joining" and treshholding the intensity and how many are matching with the targetetd mass list

Mass_List_targeted <- PseudoMALDI_targeted_extended %>% mutate(M = mass) %>% filter(M > 1000, M < 4000) %>% select(M) %>%
                      mutate(M = round(M, digits = 1)) %>% group_by(M) %>% do(head(.,1))

intensities <- c(10^4, 5*10^4, 10^5, 5*10^5, 10^6, 5*10^6, 10^7, 5*10^7, 10^8, 5*10^8, 10^9, 5*10^9, 10^10, 5*10^10) #vector with intensity thresholds
Joined_TIC <- inner_join(SNOG_Precursor_gather, SNOG_TIC, by = c("Sample", "M")) #Required data-frame with precursor intensities


###function for determination of mass bins with different intensity cutoffs
func_massbins_threshold <- function(x){
  y <- Joined_TIC %>% filter(intensity > x) %>% select(M, Sample, intensity) %>% group_by(M) %>% do(head(.,1)) ###intensity filter BEFORE group_by(M)!
  nrow(y)
}

n_mass_bins <- sapply(intensities, func_massbins_threshold) #apply function with different intensity thresholds

###function for matching with our glycoDB with different intensity cutoffs
func_massbins_threshold_matched <- function(x){
  y <- Joined_TIC %>% filter(intensity > x) %>% select(M, Sample, intensity) %>% group_by(M) %>% do(head(.,1)) ###intensity filter BEFORE group_by(M)!
  b <- inner_join(y, Mass_List_targeted, by ="M")
  nrow(b)
}

n_mass_bins_matched <- sapply(intensities, func_massbins_threshold_matched) #apply function with different intensity thresholds
mass_bins_combined <- data.frame("intensity" = intensities, "n_mass_bins" = n_mass_bins, "n_mass_bins_matched" = n_mass_bins_matched) ##combine the data

#Visualize number of mass bins and number of matched mass bins with glycoDB as a function of intensity threshold
mass_bins_combined %>% ggplot() + geom_point(aes(log(intensity), n_mass_bins), col = "red") +
                       geom_point(aes(log(intensity), n_mass_bins_matched), col = "blue") +
                       ylab("number of mass bins") + xlab("log intensity threshold") + theme_bw() 

############################################################################################################################


############################################################################################################################
#Determination of the TIC with different intensity thresholds

####function for calculating the TIC after thresholding the intensity
func_TIC_threshold <- function(x){
  y <- Joined_TIC %>% filter(intensity > x)
  sum(y$intensity)
}

TIC_raw <- sapply(intensities, func_TIC_threshold) #apply the intensity thresholds to the function
TIC_threshold <- data.frame("threshold_intensity" = intensities, "TIC" = TIC_raw, "TIC_normalized" = TIC_raw/sum(Joined_TIC$intensity))


####function for calculating the TIC after thresholding the intensity and after matching with the glycoDB
func_TIC_threshold_matched <- function(x){
  y <- Joined_TIC %>% filter(intensity > x)
  b <- inner_join(y, Mass_List_targeted, by ="M")
  sum(b$intensity)
}

TIC_matched_raw <- sapply(intensities, func_TIC_threshold_matched) #apply the intensity thresholds to the function

TIC_combined <- data.frame("intensity" = intensities, "TIC" = TIC_raw, "TIC_matched" = TIC_matched_raw) #combine the data

#Visualize TIC and TIC of matched mass bins with glycoDB as a function of intensity threshold
TIC_combined %>% ggplot() + geom_point(aes(log(intensity), TIC), col = "red") + geom_point(aes(log(intensity), TIC_matched), col = "blue") +
                 xlab("log intensity threshold") + ylab("TIC") + theme_bw() + geom_line(aes(log(intensity), TIC), alpha = 0.2) +
                 geom_line(aes(log(intensity), TIC_matched), alpha = 0.2)



#####Getting the numbers

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
matched_threshold
##############################################################################################################################################################################


library(VennDiagram)

Venn_data <- list("M_targeted" = Mass_List_targeted$M, "M_nontargeted" = SNOG_TIC_without_groups$M)

venn.diagram(Venn_data,filename = "targeted_nontargeted.tiff", fill = c("lightblue", "orange"), 
             alpha = c(0.5, 0.3), lwd =0.5, cex = 0, cat.cex = 0, label.col = 0)


Joined_TIC_all_threshold <- Joined_TIC %>% filter(intensity > 0) %>% group_by(M) %>% do(head(.,1))
Venn_data_threshold <- list("M_targeted" = targeted$M, "M_nontargeted" = Joined_TIC_all_threshold$M)

venn.diagram(Venn_data_threshold,filename = "targeted_nontargeted_threshold.tiff", fill = c("lightblue", "orange"), 
             alpha = c(0.5, 0.3), lwd =0.5, cex = 0, cat.cex = 0, label.col = 0)

############################################################################################################################################################################
###########comparison TIC targetetd vs non-targeted tissue dependent
load(file ="./RDA/PseudoMALDI_targeted_extended.rda") 

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

