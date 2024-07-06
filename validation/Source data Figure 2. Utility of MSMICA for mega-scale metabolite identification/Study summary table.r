# clear memory
rm(list = ls())

library(dplyr)
library(readr)
library(ggplot2)
library(officer)
library(rvg)

# import sample sequence list data to count the sample size
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDWB obesity/Raw data")
CHDWB_obese_mapping_hilicpos = read_delim("CHDWB_obese_mapping_hilicpos.txt")
CHDWB_obese_mapping_c18neg = read_delim("CHDWB_obese_mapping_c18neg.txt")

# Rename the column name in the sequence file: File Name -> File_Name, Sample ID -> Sample_ID
colnames(CHDWB_obese_mapping_hilicpos)[1:2] = c("File_Name", "Sample_ID")
colnames(CHDWB_obese_mapping_c18neg)[1:2] = c("File_Name", "Sample_ID")

setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDS_breast_cancer")
CHDS_mapping_hilicpos = read_delim("chds_mapping_hilicpos.txt")
CHDS_mapping_c18neg = read_delim("chds_mapping_c18neg.txt")

# Rename the column name in the sequence file: File Name -> File_Name, Sample ID -> Sample_ID
colnames(CHDS_mapping_hilicpos)[1:2] = c("File_Name", "Sample_ID")
colnames(CHDS_mapping_c18neg)[1:2] = c("File_Name", "Sample_ID")

setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDWB Urine")
CHDWB_urine_mapping_hilicpos = read_delim("chdwb_urine_mapping_hilicpos.txt")
CHDWB_urine_mapping_c18neg = read_delim("chdwb_urine_mapping_c18neg.txt")

# Rename the column name in the sequence file: File Name -> File_Name, Sample ID -> Sample_ID
colnames(CHDWB_urine_mapping_hilicpos)[1:2] = c("File_Name", "Sample_ID")
colnames(CHDWB_urine_mapping_c18neg)[1:2] = c("File_Name", "Sample_ID")

setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/RCC")
RCC_mapping_hilicpos = read_delim("rcc_mapping_hilicpos.txt")
RCC_mapping_c18neg = read_delim("rcc_mapping_c18neg.txt")

# Rename the column name in the sequence file: File Name -> File_Name, Sample ID -> Sample_ID
colnames(RCC_mapping_hilicpos)[1:2] = c("File_Name", "Sample_ID")
colnames(RCC_mapping_c18neg)[1:2] = c("File_Name", "Sample_ID")

# remove those values that are reference samples
CHDWB_obese_mapping_hilicpos = CHDWB_obese_mapping_hilicpos %>% 
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

CHDWB_obese_mapping_c18neg = CHDWB_obese_mapping_c18neg %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

CHDS_mapping_hilicpos = CHDS_mapping_hilicpos %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

CHDS_mapping_c18neg = CHDS_mapping_c18neg %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

CHDWB_urine_mapping_hilicpos = CHDWB_urine_mapping_hilicpos %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

CHDWB_urine_mapping_c18neg = CHDWB_urine_mapping_c18neg %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

RCC_mapping_hilicpos = RCC_mapping_hilicpos %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

RCC_mapping_c18neg = RCC_mapping_c18neg %>%
    filter(!(grepl("nist", Sample_ID) | grepl("q3", Sample_ID) | grepl("qstd3", Sample_ID) | grepl("q4", Sample_ID) | grepl("qstd4", Sample_ID) | grepl("SB", Sample_ID) | grepl("FB", Sample_ID)))

# count the sample size by using nrow()/3 since each sample has 3 technical replicates
CHDWB_obese_mapping_hilicpos_sample_size = nrow(CHDWB_obese_mapping_hilicpos)/3
CHDWB_obese_mapping_c18neg_sample_size = nrow(CHDWB_obese_mapping_c18neg)/3

print(CHDWB_obese_mapping_hilicpos_sample_size)
print(CHDWB_obese_mapping_c18neg_sample_size)

CHDS_mapping_hilicpos_sample_size = nrow(CHDS_mapping_hilicpos)/3
CHDS_mapping_c18neg_sample_size = nrow(CHDS_mapping_c18neg)/3

print(CHDS_mapping_hilicpos_sample_size)
print(CHDS_mapping_c18neg_sample_size)

CHDWB_urine_mapping_hilicpos_sample_size = nrow(CHDWB_urine_mapping_hilicpos)/3
CHDWB_urine_mapping_c18neg_sample_size = nrow(CHDWB_urine_mapping_c18neg)/3

print(CHDWB_urine_mapping_hilicpos_sample_size)
print(CHDWB_urine_mapping_c18neg_sample_size)

RCC_mapping_hilicpos_sample_size = nrow(RCC_mapping_hilicpos)/3
RCC_mapping_c18neg_sample_size = nrow(RCC_mapping_c18neg)/3

print(RCC_mapping_hilicpos_sample_size)
print(RCC_mapping_c18neg_sample_size)


# import MSMICA data for each column mode 
## CHDWB obesity
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDWB obesity/MSMICA/triplicate/CHDWB_obesity_MSMICA_validation")
## HILIC pos
CHDWB_hilicpos_MSMICA_identified_metabolites = read_csv("CHDWB_obesity_hilicpos_MSMICA_identified_metabolites.csv")
## C18 neg
CHDWB_c18neg_MSMICA_identified_metabolites = read_csv("CHDWB_obesity_c18neg_MSMICA_identified_metabolites.csv")

## CHDS breast cancer
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDS_breast_cancer/MSMICA/triplicate/CHDS_MSMICA_validation")
## HILIC pos
CHDS_hilicpos_MSMICA_identified_metabolites = read_csv("CHDS_hilicpos_MSMICA_identified_metabolites.csv")
## C18 neg
CHDS_c18neg_MSMICA_identified_metabolites = read_csv("CHDS_c18neg_MSMICA_identified_metabolites.csv")

## CHDWB urine
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDWB Urine/triplicate/CHDWB_urine_MSMICA_validation")
## HILIC pos
CHDWB_urine_hilicpos_MSMICA_identified_metabolites = read_csv("CHDWB_urine_hilicpos_MSMICA_identified_metabolites.csv")
## C18 neg
CHDWB_urine_c18neg_MSMICA_identified_metabolites = read_csv("CHDWB_urine_c18neg_MSMICA_identified_metabolites.csv")

## SLAM
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/SLAM/Annotation - Raw/batch16/triplicate/SLAM_MSMICA_validation")
## HILIC pos
SLAM_hilicpos_MSMICA_identified_metabolites = read_csv("SLAM_hilicpos_MSMICA_identified_metabolites.csv")
## C18 neg
SLAM_c18neg_MSMICA_identified_metabolites = read_csv("SLAM_c18neg_MSMICA_identified_metabolites.csv")

## RCC
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/RCC/triplicate/RCC_MSMICA_validation")
## HILIC pos
RCC_hilicpos_MSMICA_identified_metabolites = read_csv("RCC_target_hilicpos_MSMICA_identified_metabolites.csv")
## C18 neg
RCC_c18neg_MSMICA_identified_metabolites = read_csv("RCC_target_c18neg_MSMICA_identified_metabolites.csv")

# count the number of unique identified metabolites by MSMICA using KEGGID column in each study, hilicpos, c18neg, and total
## CHDWB obesity
CHDWB_hilicpos_MSMICA_identified_metabolites_unique = length(unique(CHDWB_hilicpos_MSMICA_identified_metabolites$KEGGID)) 
CHDWB_c18neg_MSMICA_identified_metabolites_unique = length(unique(CHDWB_c18neg_MSMICA_identified_metabolites$KEGGID))
CHDWB_total_MSMICA_identified_metabolites_unique = length(unique(c(CHDWB_hilicpos_MSMICA_identified_metabolites$KEGGID, CHDWB_c18neg_MSMICA_identified_metabolites$KEGGID)))

print(CHDWB_hilicpos_MSMICA_identified_metabolites_unique)
print(CHDWB_c18neg_MSMICA_identified_metabolites_unique)
print(CHDWB_total_MSMICA_identified_metabolites_unique)

## CHDS breast cancer
CHDS_hilicpos_MSMICA_identified_metabolites_unique = length(unique(CHDS_hilicpos_MSMICA_identified_metabolites$KEGGID))
CHDS_c18neg_MSMICA_identified_metabolites_unique = length(unique(CHDS_c18neg_MSMICA_identified_metabolites$KEGGID))
CHDS_total_MSMICA_identified_metabolites_unique = length(unique(c(CHDS_hilicpos_MSMICA_identified_metabolites$KEGGID, CHDS_c18neg_MSMICA_identified_metabolites$KEGGID)))

print(CHDS_hilicpos_MSMICA_identified_metabolites_unique)
print(CHDS_c18neg_MSMICA_identified_metabolites_unique)
print(CHDS_total_MSMICA_identified_metabolites_unique)

## CHDWB urine
CHDWB_urine_hilicpos_MSMICA_identified_metabolites_unique = length(unique(CHDWB_urine_hilicpos_MSMICA_identified_metabolites$KEGGID))
CHDWB_urine_c18neg_MSMICA_identified_metabolites_unique = length(unique(CHDWB_urine_c18neg_MSMICA_identified_metabolites$KEGGID))
CHDWB_urine_total_MSMICA_identified_metabolites_unique = length(unique(c(CHDWB_urine_hilicpos_MSMICA_identified_metabolites$KEGGID, CHDWB_urine_c18neg_MSMICA_identified_metabolites$KEGGID)))

print(CHDWB_urine_hilicpos_MSMICA_identified_metabolites_unique)
print(CHDWB_urine_c18neg_MSMICA_identified_metabolites_unique)
print(CHDWB_urine_total_MSMICA_identified_metabolites_unique)

## SLAM
SLAM_hilicpos_MSMICA_identified_metabolites_unique = length(unique(SLAM_hilicpos_MSMICA_identified_metabolites$KEGGID))
SLAM_c18neg_MSMICA_identified_metabolites_unique = length(unique(SLAM_c18neg_MSMICA_identified_metabolites$KEGGID))
SLAM_total_MSMICA_identified_metabolites_unique = length(unique(c(SLAM_hilicpos_MSMICA_identified_metabolites$KEGGID, SLAM_c18neg_MSMICA_identified_metabolites$KEGGID)))

print(SLAM_hilicpos_MSMICA_identified_metabolites_unique)
print(SLAM_c18neg_MSMICA_identified_metabolites_unique)
print(SLAM_total_MSMICA_identified_metabolites_unique)

## RCC
RCC_hilicpos_MSMICA_identified_metabolites_unique = length(unique(RCC_hilicpos_MSMICA_identified_metabolites$KEGGID))
RCC_c18neg_MSMICA_identified_metabolites_unique = length(unique(RCC_c18neg_MSMICA_identified_metabolites$KEGGID))
RCC_total_MSMICA_identified_metabolites_unique = length(unique(c(RCC_hilicpos_MSMICA_identified_metabolites$KEGGID, RCC_c18neg_MSMICA_identified_metabolites$KEGGID)))

print(RCC_hilicpos_MSMICA_identified_metabolites_unique)
print(RCC_c18neg_MSMICA_identified_metabolites_unique)
print(RCC_total_MSMICA_identified_metabolites_unique)

# import all validation results for all studies
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDWB obesity/MSMICA/triplicate/CHDWB_obesity_MSMICA_validation")
CHDWB_hilicpos_validation = read_csv("validation_target_hilicpos_MSMICA_identified_metabolites.csv")
CHDWB_c18neg_validation = read_csv("validation_target_c18neg_MSMICA_identified_metabolites.csv")
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/SLAM/Annotation - Raw/batch16/triplicate/SLAM_MSMICA_validation")
SLAM_hilicpos_validation = read_csv("validation_target_hilicpos_MSMICA_identified_metabolites.csv")
SLAM_c18neg_validation = read_csv("validation_target_c18neg_MSMICA_identified_metabolites.csv")
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDS_breast_cancer/MSMICA/triplicate/CHDS_MSMICA_validation")
CHDS_hilicpos_validation = read_csv("validation_target_hilicpos_MSMICA_identified_metabolites.csv")
CHDS_c18neg_validation = read_csv("validation_target_c18neg_MSMICA_identified_metabolites.csv")
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/CHDWB Urine/triplicate/CHDWB_urine_MSMICA_validation")
CHDWB_urine_hilicpos_validation = read_csv("validation_target_hilicpos_MSMICA_identified_metabolites.csv")
CHDWB_urine_c18neg_validation = read_csv("validation_target_c18neg_MSMICA_identified_metabolites.csv")
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/RCC/triplicate/RCC_MSMICA_validation")
RCC_hilicpos_validation = read_csv("validation_target_hilicpos_MSMICA_identified_metabolites.csv")
RCC_c18neg_validation = read_csv("validation_target_c18neg_MSMICA_identified_metabolites.csv")

nrow(CHDWB_hilicpos_validation)
nrow(CHDWB_c18neg_validation)
nrow(SLAM_hilicpos_validation)
nrow(SLAM_c18neg_validation)
nrow(CHDS_hilicpos_validation)
nrow(CHDS_c18neg_validation)
nrow(CHDWB_urine_hilicpos_validation)
nrow(CHDWB_urine_c18neg_validation)
nrow(RCC_hilicpos_validation)
nrow(RCC_c18neg_validation)


# extract the Fraction_correct and False_discovery_rate values from the first row of each validation result
CHDWB_hilicpos_validation_2 = tibble(
    Study = "CHDWB obesity HILIC positive",
    ## count the row number where the match column is "match" to get the true positive
    True_positive = nrow(CHDWB_hilicpos_validation[CHDWB_hilicpos_validation$match == "match", ]),
    ## count the row number where the match column is not "match" to get the false positive
    False_positive = nrow(CHDWB_hilicpos_validation[CHDWB_hilicpos_validation$match != "match", ]),
    Total = nrow(CHDWB_hilicpos_validation),
    Fraction_correct = CHDWB_hilicpos_validation$Fraction_correct[1],
    False_discovery_rate = CHDWB_hilicpos_validation$False_discovery_rate[1]
)

CHDWB_c18neg_validation_2 = tibble(
    Study = "CHDWB obesity C18 negative",
    True_positive = nrow(CHDWB_c18neg_validation[CHDWB_c18neg_validation$match == "match", ]),
    False_positive = nrow(CHDWB_c18neg_validation[CHDWB_c18neg_validation$match != "match", ]),
    Total = nrow(CHDWB_c18neg_validation),
    Fraction_correct = CHDWB_c18neg_validation$Fraction_correct[1],
    False_discovery_rate = CHDWB_c18neg_validation$False_discovery_rate[1]
)

SLAM_hilicpos_validation_2 = tibble(
    Study = "SLAM HILIC positive",
    True_positive = nrow(SLAM_hilicpos_validation[SLAM_hilicpos_validation$match == "match", ]),
    False_positive = nrow(SLAM_hilicpos_validation[SLAM_hilicpos_validation$match != "match", ]),
    Total = nrow(SLAM_hilicpos_validation),
    Fraction_correct = SLAM_hilicpos_validation$Fraction_correct[1],
    False_discovery_rate = SLAM_hilicpos_validation$False_discovery_rate[1]
)

SLAM_c18neg_validation_2 = tibble(
    Study = "SLAM C18 negative",
    True_positive = nrow(SLAM_c18neg_validation[SLAM_c18neg_validation$match == "match", ]),
    False_positive = nrow(SLAM_c18neg_validation[SLAM_c18neg_validation$match != "match", ]),
    Total = nrow(SLAM_c18neg_validation),
    Fraction_correct = SLAM_c18neg_validation$Fraction_correct[1],
    False_discovery_rate = SLAM_c18neg_validation$False_discovery_rate[1]
)

CHDS_hilicpos_validation_2 = tibble(
    Study = "CHDS breast cancer HILIC positive",
    True_positive = nrow(CHDS_hilicpos_validation[CHDS_hilicpos_validation$match == "match", ]),
    False_positive = nrow(CHDS_hilicpos_validation[CHDS_hilicpos_validation$match != "match", ]),
    Total = nrow(CHDS_hilicpos_validation),
    Fraction_correct = CHDS_hilicpos_validation$Fraction_correct[1],
    False_discovery_rate = CHDS_hilicpos_validation$False_discovery_rate[1]
)

CHDS_c18neg_validation_2 = tibble(
    Study = "CHDS breast cancer C18 negative",
    True_positive = nrow(CHDS_c18neg_validation[CHDS_c18neg_validation$match == "match", ]),
    False_positive = nrow(CHDS_c18neg_validation[CHDS_c18neg_validation$match != "match", ]),
    Total = nrow(CHDS_c18neg_validation),
    Fraction_correct = CHDS_c18neg_validation$Fraction_correct[1],
    False_discovery_rate = CHDS_c18neg_validation$False_discovery_rate[1]
)

CHDWB_urine_hilicpos_validation_2 = tibble(
    Study = "CHDWB urine HILIC positive",
    True_positive = nrow(CHDWB_urine_hilicpos_validation[CHDWB_urine_hilicpos_validation$match == "match", ]),
    False_positive = nrow(CHDWB_urine_hilicpos_validation[CHDWB_urine_hilicpos_validation$match != "match", ]),
    Total = nrow(CHDWB_urine_hilicpos_validation),
    Fraction_correct = CHDWB_urine_hilicpos_validation$Fraction_correct[1],
    False_discovery_rate = CHDWB_urine_hilicpos_validation$False_discovery_rate[1]
)

CHDWB_urine_c18neg_validation_2 = tibble(
    Study = "CHDWB urine C18 negative",
    True_positive = nrow(CHDWB_urine_c18neg_validation[CHDWB_urine_c18neg_validation$match == "match", ]),
    False_positive = nrow(CHDWB_urine_c18neg_validation[CHDWB_urine_c18neg_validation$match != "match", ]),
    Total = nrow(CHDWB_urine_c18neg_validation),
    Fraction_correct = CHDWB_urine_c18neg_validation$Fraction_correct[1],
    False_discovery_rate = CHDWB_urine_c18neg_validation$False_discovery_rate[1]
)

RCC_hilicpos_validation_2 = tibble(
    Study = "RCC HILIC positive",
    True_positive = nrow(RCC_hilicpos_validation[RCC_hilicpos_validation$match == "match", ]),
    False_positive = nrow(RCC_hilicpos_validation[RCC_hilicpos_validation$match != "match", ]),
    Total = nrow(RCC_hilicpos_validation),
    Fraction_correct = RCC_hilicpos_validation$Fraction_correct[1],
    False_discovery_rate = RCC_hilicpos_validation$False_discovery_rate[1]
)

RCC_c18neg_validation_2 = tibble(
    Study = "RCC C18 negative",
    True_positive = nrow(RCC_c18neg_validation[RCC_c18neg_validation$match == "match", ]),
    False_positive = nrow(RCC_c18neg_validation[RCC_c18neg_validation$match != "match", ]),
    Total = nrow(RCC_c18neg_validation),
    Fraction_correct = RCC_c18neg_validation$Fraction_correct[1],
    False_discovery_rate = RCC_c18neg_validation$False_discovery_rate[1]
)

# combine all validation results
MSMICA_validation_result = bind_rows(CHDWB_hilicpos_validation_2, CHDWB_c18neg_validation_2, SLAM_hilicpos_validation_2, SLAM_c18neg_validation_2, CHDS_hilicpos_validation_2, CHDS_c18neg_validation_2, CHDWB_urine_hilicpos_validation_2, CHDWB_urine_c18neg_validation_2, RCC_hilicpos_validation_2, RCC_c18neg_validation_2)

# calculate the average fraction correct and false discovery rate
MSMICA_validation_result_summary = MSMICA_validation_result %>%
    summarise(
        True_positive = sum(True_positive),
        False_positive = sum(False_positive),
        Total = sum(Total),
        Fraction_correct = True_positive/Total,
        False_discovery_rate = False_positive/Total
    ) %>%
    mutate(Study = "Average") %>%
    select(Study, True_positive, False_positive, Total, Fraction_correct, False_discovery_rate)

# add the average fraction correct and false discovery rate to the validation result
MSMICA_validation_result = bind_rows(MSMICA_validation_result, MSMICA_validation_result_summary)

# factorize the Study column
MSMICA_validation_result$Study = factor(MSMICA_validation_result$Study, levels = c("CHDWB obesity HILIC positive", "CHDWB obesity C18 negative", "CHDS breast cancer HILIC positive", "CHDS breast cancer C18 negative", "CHDWB urine HILIC positive", "CHDWB urine C18 negative", "SLAM HILIC positive", "SLAM C18 negative", "RCC HILIC positive", "RCC C18 negative", "Average"))

# save the summary result
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/Publication/Abstract/MSMICA_summary")
write_csv(MSMICA_validation_result, "MSMICA_validation_result_summary.csv")

# use ggplot to plots the validation results: x-axis = Study, y-axis = Fraction_correct, color = Study using geom_col
MSMICA_validation_fraction_correct_plot = ggplot(MSMICA_validation_result, aes(x = Study, y = Fraction_correct, fill = Study)) +
    geom_col() +
    labs(y = "Fraction correct") +
    # minimize the theme
    theme_minimal() +
    # hide the x axis label and title
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),  # Increase font size for y-axis title
        axis.text.y = element_text(size = 12),   # Increase font size for y-axis text
        legend.title = element_text(size = 12),  # Increase font size for legend title
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        ) +
    # make y axis from 0.7 to 1 with 0.05 interval
    coord_cartesian(ylim=c(0.7, 1)) +
    scale_y_continuous(breaks = seq(0.7, 1, 0.05))

# use ggplot to plots the validation results: x-axis = Study, y-axis = False_discovery_rate, color = Study using geom_col
MSMICA_validation_false_discovery_rate_plot = ggplot(MSMICA_validation_result, aes(x = Study, y = False_discovery_rate, fill = Study)) +
    geom_col() +
    labs(y = "False discovery rate") +
    # minimize the theme
    theme_minimal() +
    # hide the x axis label and title
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),  # Increase font size for y-axis title
        axis.text.y = element_text(size = 12),   # Increase font size for y-axis text
        legend.title = element_text(size = 12),  # Increase font size for legend title
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        ) +
    # make y axis from 0 to 0.1 with 0.02 interval
    coord_cartesian(ylim=c(0, 0.1)) +
    scale_y_continuous(breaks = seq(0, 0.1, 0.02))


######################
# Convert ggplot objects to dml objects
MSMICA_validation_fraction_correct_plot_vec <- dml(ggobj = MSMICA_validation_fraction_correct_plot, fonts = list(sans = "Times New Roman"))
MSMICA_validation_false_discovery_rate_plot_vec <- dml(ggobj = MSMICA_validation_false_discovery_rate_plot, fonts = list(sans = "Times New Roman"))

# Create a new empty pptx file
doc <- read_pptx()

# Add the slide 
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, MSMICA_validation_fraction_correct_plot_vec, ph_location_fullsize())

doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, MSMICA_validation_false_discovery_rate_plot_vec, ph_location_fullsize())

print(doc, target = "MSMICA_validation_fraction_correct_FDR_plot.pptx")