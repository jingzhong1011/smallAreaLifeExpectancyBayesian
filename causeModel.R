library(INLA)
library(dplyr)
library(data.table)
library(ggplot2)

setwd("/Users/jing-zhongwang/Library/CloudStorage/GoogleDrive-clockwcc@gmail.com/我的雲端硬碟/NTUEPM/Master Thesis/Paper Writing/Data")

DeathCount_cause <- fread("DeathCount_cause.csv")

Indi <- read.csv("MOIMOH_ZH.csv", fileEncoding = "BIG5")

DC_cause <- DeathCount_cause %>% left_join(Indi, by = "MOI") %>% 
  group_by(YEAR, SEX, AGP, CAUSE, Indi) %>% 
  summarise(COUNT = sum(COUNT),
            POP = sum(POP)) %>% 
  ungroup() %>% 
  mutate(PROP = COUNT / POP)

DC_cause_male <- DC_cause %>% filter(SEX == 1)
DC_cause_female <- DC_cause %>% filter(SEX == 2)
unique(DC_cause$CAUSE)


#######
causeModeling <- function(data){
  results <- list()
  sex_idx = c(1: 2)
  cause_idx <- unique(data$CAUSE)
  
  for(i in sex_idx){
    for(j in cause_idx){
      cat("Processing: SEX =", i, "CAUSE =", j, "\n")
      dt <- data %>% filter(SEX == i & CAUSE == j) 
      dt2 <- dt %>% 
        mutate(AGP2 = AGP,
               AGP3 = AGP,
               YEAR2 = YEAR,
               YEAR3 = YEAR,
               Indi2 = Indi,
               Indi3 = Indi)

      if(sum(dt2$COUNT) == 0){
        dt_fitted <- dt %>% mutate(PROP_SMOOTHED = 0)
        results[[paste(i, j, sep = "_")]] <- dt_fitted
        next  
      }
      fit <- inla(COUNT ~ 
                    f(AGP, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(YEAR, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(Indi, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(YEAR2, AGP2, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(Indi2, AGP3, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(YEAR3, Indi3, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))),
           family = "poisson", data = dt2, E = POP, control.predictor = list(link = 1))
      
      dt_fitted <- dt %>% mutate(PROP_SMOOTHED = fit$summary.fitted.values$mean)
      results[[paste(i, j, sep = "_")]] <- dt_fitted
    }
  }
 return(results) 
}
output <- causeModeling(DC_cause)


output2 <- bind_rows(output)


write.csv(output2, "DCcauseSmoothed.csv", row.names = F)

##########################

stand_pop <- data.frame(AGP = c(0, 1, seq(5, 85, by = 5)),
                        Percent = c(17917, 70652, 86870, 85970, 84670,
                                    82171, 79272, 76073, 71475, 65877,
                                    60379, 53681, 45484, 37187, 29590,
                                    22092, 15195, 9097, 6348)) %>% mutate(Percent = Percent / 1e06)



output2_stand <- output2 %>% left_join(stand_pop, by = "AGP") %>% 
  mutate(PROP_SMOOTHED_STAND = PROP_SMOOTHED * Percent * 100000,
         PROP_STAND = PROP * Percent * 100000) %>% 
  group_by(YEAR, SEX, Indi, CAUSE) %>% 
  summarise(PROP_SMOOTHED_STAND = sum(PROP_SMOOTHED_STAND),
            PROP_STAND = sum(PROP_STAND)) %>% ungroup() %>% 
  mutate(Indi_eng = factor(case_when(Indi == 0 ~ "Non-Indi.",
                              Indi == 1 ~ "Mountain Indi.",
                              Indi == 2 ~ "Plain Indi."), levels = c("Non-Indi.", "Plain Indi.", "Mountain Indi.")))


ggplot(data = output2_stand %>% filter(SEX == 1)) +
  geom_line(aes(x = YEAR, y = PROP_SMOOTHED_STAND, group = Indi_eng, colour = Indi_eng)) + 
  geom_point(aes(x = YEAR, y = PROP_STAND, group = Indi_eng, colour = Indi_eng), pch = 1) +
  theme_bw() +
  facet_wrap(CAUSE ~ ., scales = "free") +
  xlab("Year") +
  ylab("Age-standardized Mortality Rate (per 100,000 population)") +
  ggtitle("Male") + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = output2_stand %>% filter(SEX == 2)) +
  geom_line(aes(x = YEAR, y = PROP_SMOOTHED_STAND, group = Indi_eng, colour = Indi_eng)) + 
  geom_point(aes(x = YEAR, y = PROP_STAND, group = Indi_eng, colour = Indi_eng), pch = 1) +
  theme_bw() +
  facet_wrap(CAUSE ~ ., scales = "free") +
  xlab("Year") +
  ylab("Age-standardized Mortality Rate (per 100,000 population)") +
  ggtitle("Female") + 
  theme(plot.title = element_text(hjust = 0.5))
