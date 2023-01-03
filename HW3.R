library(tidyverse)
library(magrittr)
library(haven)
library(gmodels)
library(psych)
library(xtable)


#----------  Set up theme
#===============================================================================
mytheme <- function(...) {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14,color = "grey10",  face = "bold", hjust = 0.5),
      plot.subtitle = element_text(face = "italic", color = "gray10", size = 14),
      plot.caption = element_text(face = "italic", size = 14, color = "gray10"),
      axis.line = element_line(linetype = "solid"),
      axis.text.x = element_text(color = "gray10", size = 10),
      axis.text.y = element_text(color = "gray10", size = 10),
      # axis.ticks = element_blank(),
      axis.title.x = element_text(color = "gray10", size = 14),
      axis.title.y = element_text(color = "gray10", size = 14),
      panel.grid.minor = element_blank(),
      # panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal",
      legend.position = c(0.5, 0.9),
      legend.background = element_rect(fill = NA, color = NA),
      legend.text = element_text(size = 12),
      legend.key.width = unit(2, "line"),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = NA, color = NA)
    )
}

df_org <- read_sas("hearing500lr.sas7bdat")

# Round time as requested
df_org %<>% mutate(timeround = round(TIME))


df_org %<>% group_by(id, side, timeround) %>% 
    summarise(y = mean(y),
              age = mean(age)) %>% ungroup() %>%
    mutate(ycat = ifelse(y <= 15, 1, ifelse(y <= 25, 2, 3)),
           side = case_when(side == "left" ~ "Left",
                            side == "right" ~ "Right"))
df_org

# Extract age for merging
age_uniqe <- df_org %>% group_by(id) %>% summarise(age = mean(age)) %>% ungroup()


# Create frame for merging
df_frame <- expand.grid(id = unique(df_org$id), side = c("Left", "Right"), 
                        timeround = c(0:22)) %>% as_tibble()

# Merge data to frame
# 546*23*2 = 25116

df <- df_frame %>% left_join(df_org, by = c("id", "side", "timeround")) %>%
    arrange(id, timeround, side) %>% select(-age) %>% left_join(age_uniqe, by = "id")

df

# 25116 - 4407 = 20709
sapply(df, function(x){sum(is.na(x))})


# df %<>% group_by(id, side) %>%
#     mutate(ylag = lag(y, 1)) %>% ungroup()

# saveRDS(df, "hearing_timeround.RDS")
# write.csv(df, "hearing_timeround.csv", row.names = F, na = "")

# df %>% select(id, side, timeround, y, ylag) %>% View()




df_wide <- pivot_wider(df %>% select(-ycat), names_from = timeround, values_from = y) 
names(df_wide) <- c("id", "side", "age", paste0("t", 0:22))

# saveRDS(df_wide, "hearing_timeround_wide.RDS")
# write.csv(df_wide, "hearing_timeround_wide.csv", row.names = F, na = "")


sapply(df_wide, function(x){sum(is.na(x))})



hist(df$y)

set.seed(12345)
df_org %>% filter(id %in% sample(unique(df_org$id), 30)) %>%
    ggplot(aes(x = as.factor(timeround), y = y, group = id)) +
    geom_point() +
    geom_line() +
    facet_wrap(~side, ncol = 1) +
    mytheme()


df %>%
    ggplot(aes(x = as.factor(timeround), y = y, group = id)) +
    geom_point() +
    geom_line() +
    facet_wrap(~side, ncol = 1) +
    mytheme()


df_org %>% group_by(side, timeround) %>%
    summarise(y_mean = mean(y),
              y_sd= sd(y),
              n = n(),
              se = y_sd/sqrt(n)) %>%
    ggplot(aes(x = as.factor(timeround), y = y_mean, fill = side, color = side, group = side)) +
    geom_point(size = 3, shape = 21, show.legend = F) +
    geom_line(linewidth = 1, show.legend = F) +
    geom_errorbar(aes(ymin = y_mean - se, ymax = y_mean + se),
                  colour="black", width=0.2, linewidth = 0.5) +
    facet_wrap(~side, ncol = 1) +
    mytheme()



x <- df_org %>% group_by(side, timeround) %>%
    summarise(n = n())

x <- pivot_wider(x, names_from = timeround, values_from = n)

# xtable(x)
x
apply(x[, 2:24], 1, sum)


#---------- Missing value pattern
#===============================================================================
df_wide$mdpattern <- NA

for (i in 1:nrow(df_wide)) {
    pattern <- NA
    if (!is.na(df_wide[i, 4])) {
        pattern <- "1"
    } else {
        pattern <- "0"
    }
    for (j in 5:26) {
        if (!is.na(df_wide[i, j])) {
            pattern <- paste0(pattern, "1")
        } else {
            pattern <- paste0(pattern, "0")
        }
    }
    df_wide[i, 27] <- pattern
}


df_wide %>% select(mdpattern) %>% View()

monotone <- c(
    "00000000000000000000000",
    "10000000000000000000000",
    "11000000000000000000000",
    "11100000000000000000000",
    "11110000000000000000000",
    "11111000000000000000000",
    "11111100000000000000000",
    "11111110000000000000000",
    "11111111000000000000000",
    "11111111100000000000000",
    "11111111110000000000000",
    "11111111111000000000000",
    "11111111111100000000000",
    "11111111111110000000000",
    "11111111111111000000000",
    "11111111111111100000000",
    "11111111111111110000000",
    "11111111111111111000000",
    "11111111111111111100000",
    "11111111111111111110000",
    "11111111111111111111000",
    "11111111111111111111100",
    "11111111111111111111110",
    "11111111111111111111111"
)


mi_matrix_l <- cbind(monotone, n = NA) %>% as.matrix()

for (i in 1:24) {
    mi_matrix_l[i, 2] <- df_wide %>% filter(mdpattern %in% mi_matrix_l[i, 1] & side == "Left") %>%
        nrow()
}

mi_matrix_r <- cbind(monotone, n = NA) %>% as.matrix()
for (i in 1:24) {
    mi_matrix_r[i, 2] <- df_wide %>% filter(mdpattern %in% mi_matrix_r[i, 1] & side == "Right") %>%
        nrow()
}


# Percentage of monotone missingness
df_wide %>% filter(mdpattern %in% monotone & side == "Left") %>%
    nrow()/546

df_wide %>% filter(mdpattern %in% monotone & side == "Right") %>%
    nrow()/546


# Using package
library(finalfit)

df_wide_l <- df_wide %>% filter(side == "Left") %>% arrange(id)
df_wide_r <- df_wide %>% filter(side == "Right") %>% arrange(id)

md_l <- missing_plot(df_wide_l[, 4:26], title = "Left ear") + 
    theme(plot.title = element_text(size = 16,color = "grey10",  face = "bold", hjust = 0.5))
md_r <- missing_plot(df_wide_r[, 4:26], title = "Right ear") +
    theme(plot.title = element_text(size = 16,color = "grey10",  face = "bold", hjust = 0.5))

gridExtra::grid.arrange(md_l, md_r, ncol = 2)


library(naniar)
mcar_test(df_wide_l[, 4:26])
mcar_test(df_wide_r[, 4:26])




#---------- Q4
#===============================================================================
# Reshape data
df_wide_im <- read_sas("hearing_wide_im.sas7bdat")

df_wide_im %<>% gather(-c("_Imputation_", "id", "side", "age", "side2"), 
                       key = "time", value = "y")

df_wide_im %<>% mutate(time = parse_number(time))


glimpse(df_wide_im)


table(df_wide_im$time)

# write.csv(df_wide_im, "df_long_im.csv", row.names = F, na = "")









