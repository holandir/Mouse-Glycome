###Impact of different SNOG scores on number of MSMS spectra and unique mass bins
library(tidyverse)

#load files
load(file = "./RDA/PseudoMALDI_targeted_extended.rda")
load(file= "./RDA/SNOG_complete_diagnostic.rda")

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