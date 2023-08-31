##########################################################################################################################################
library(tidyverse)

load(file = "./RDA/PseudoMALDI_targeted_extended.rda")
load(file= "./RDA/SNOG_complete_diagnostic.rda")
load(file= "./RDA/SNOG_Precursor_gather.rda")
load(file= "./RDA/Joined_TIC.rda")


#######How many mass bins are there after "joining" and treshholding the intensity and how many are matching with the targetetd mass list

###Mass-list from glycoDB
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
