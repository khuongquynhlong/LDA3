library(tidyverse)
library(magrittr)
library(haven)

# Prepare data for proc mianalyze
Imputation <- rep(c(1:10), each = 2)
Parameter <- rep(c("Prm6", "Prm7"), 10)
Effect <- rep(c("Prm6", "Prm7"), 10)
gmpinfoint <- tibble(`_Imputation_` = Imputation, Parameter, Effect)


df_gmcovb <- read_sas("gmcovb.sas7bdat") %>%
    rename("Prm6" = Intercept1,
           "Prm7" = Intercept2) %>%
    mutate(RowName = case_when(RowName == "Intercept1" ~ "Prm6",
                               RowName == "Intercept2" ~ "Prm7",
                               TRUE ~ RowName))

df_gmparms <- read_sas("gmparms.sas7bdat") %>%
    mutate(Parm = case_when(Parm == "Intercept1" ~ "Prm6",
                            Parm == "Intercept2" ~ "Prm7",
                            TRUE ~ Parm))

df_gmpinfo <- read_sas("gmpinfo.sas7bdat") %>% 
    add_row(gmpinfoint) %>% arrange(`_Imputation_`, Parameter)


write.csv(df_gmcovb, "gmcovb.csv", row.names = F)
write.csv(df_gmparms, "gmparms.csv", row.names = F)
write.csv(df_gmpinfo, "gmpinfo.csv", row.names = F)



