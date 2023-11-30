# Import libraries (install first if required)
library(tidyverse)
library(DescTools)
library(dplyr)

##############################
# USER DEFINED VARIABLES
# Set the working directory
setwd("path/to/working/directory/")
#Define the dictionary (list of the wells used)
dictionary <- read.csv("Dictionary.csv")
# Define the csv containing the parsed microplate data
df <- read_csv("reformatted.csv")
# Minimum Vp value for well to be considered to contain phage
Vp_min <- 0.5
###############################

# Create well factor
df$well <- factor(df$well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
df <- df %>% mutate(row=gsub('^([A-H])(.*)', '\\1', well), column=factor(gsub('^([A-H])(.*)', '\\2', well), levels=1:12))

# Include dictionary data and remove empty wells
df <- left_join(df, dictionary, by = "well")
df <- df %>% filter(isolation != "nothing")

# Normalise growth curves based on the blanks
blk_df <- df %>% filter(isolation == "blk") %>% group_by(time_hours) %>% summarise(mean_blk_od = mean(mean_od))
df <- df %>% left_join(blk_df) %>% mutate(norm_od = mean_od - mean_blk_od)

# Calculate virulence index
control <- df %>% filter(isolation == "host") %>% mutate(cntrlAUC = AUC(x = time_hours, y = mean_od, method = "trapezoid")) %>% select(well, cntrlAUC) %>% unique()
df_AUC <- df %>% group_by(well) %>% mutate(phageAUC = AUC(x = time_hours, y = mean_od, method = "trapezoid"))
df_AUC <- df_AUC %>% mutate(cntrlAUC = as.numeric(control$cntrlAUC)) %>% mutate(Vp = (1 - (phageAUC / cntrlAUC)))

df <- df %>% left_join(df_AUC)
df <- df %>% mutate(Vp = round(Vp, digits = 3))

df_automate <- df %>% select(well, time_hours, isolation, norm_od, Vp) %>% unique()

# Explicitly set the order of well levels
well_order <- LETTERS[1:8] %>% rep(each = 12) %>% paste0(1:12)
df_automate$well <- factor(df_automate$well, levels = well_order)

# Create borders to highlight phage wells in plot
df_automate <- df_automate %>%
  mutate(border = ifelse(Vp > Vp_min & isolation != "blk", "Phage", "No Phage")) %>%
  mutate(border = as.character(border)) %>%
  mutate(border = factor(border, levels = c("Phage", "No Phage")))

df_final <- df_automate %>% select(well, border, isolation) %>% unique() %>% mutate(time_hours = 0) %>% mutate(norm_od = 0)

# Create plot
Vp_plot <- ggplot(data = df_automate, aes(x = time_hours, y = norm_od)) +
  geom_rect(data = df_final, aes(fill=border), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.3) +
  geom_line(size = 2) +
  facet_wrap(~well + isolation, ncol = 12) +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600', limits = c(-0.3,0.9), breaks = seq(0, 0.8, by = 0.2)) +
  scale_fill_manual(values = c("#ffad00", "white")) +
  labs(fill = " ") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
Vp_plot

# Output selected phages
picked <- df_automate %>% filter(border == "Phage" & isolation != "blk" & isolation != "host") %>% select(well, Vp) %>% unique() %>% as.list()

write.csv(picked, "test_list_phages.csv", row.names = FALSE)
print("CSV file 'test_list_phages.csv' has been successfully created.")
