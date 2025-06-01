setwd("C:/Canteen/Western/2020 Summer/COVID19/Plot/Infection rate")
library(tidyverse)

## Import and attach the data for Ontario
library(RCurl)
myurl <- "https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv"
mydat <- read.csv(myurl, header=TRUE)
mydat = mydat[mydat$Client_Gender == 'MALE' | mydat$Client_Gender == 'FEMALE', ]
ON_case = mydat[, which(colnames(mydat) %in% c('Client_Gender', 'Age_Group'))]

## Remove unused levels
ON_case = ON_case[ON_case$Age_Group != 'UNKNOWN', ]
ON_case = droplevels(ON_case)

## Read Ontario population data
ON_pop = read.csv('ON_pop.csv', stringsAsFactors = T)
ON_pop$Age = ON_pop$ï..Age
ON_pop = droplevels(ON_pop)


## Reassign levels
ON_case <- ON_case[!ON_case$Age_Group == '', ]
ON_case = droplevels(ON_case)
levels(ON_case$Age_Group) = c('0 - 19 ', '20 - 29 ', '30 - 39 ',
                     '40 - 49 ', '50 - 59 ', '60 - 69 ', 
                     '70 - 79 ', '80 - 89 ', '90 - 99 ')


levels(ON_pop$Age) = c(rep('0 - 19 ', 3),
                    rep('20 - 29 ', 2), rep('30 - 39 ', 2), 
                    rep('40 - 49 ', 2), '0 - 19 ', 
                    rep('50 - 59 ', 2), 
                    rep('60 - 69 ', 2), rep('70 - 79 ', 2), 
                    rep('80 - 89 ', 2), rep('90 - 99 ', 2))

# Remove commas from the number
library(tidyverse)
ON_pop$Count = as.numeric(gsub(",","",ON_pop$Count))
ON_pop_count = ON_pop %>% 
  group_by(Age, Gender) %>% 
  summarise(Total = sum(Count))

ON_case_count = as.data.frame(table(ON_case))
levels(ON_case_count$Age_Group) = c('0 - 19 ', '20 - 29 ', '30 - 39 ',
                                    '40 - 49 ', '50 - 59 ', '60 - 69 ', 
                                    '70 - 79 ', '80 - 89 ', '90 - 99 ')
colnames(ON_case_count) = c('Age', 'Gender', 'Freq')
join_count = merge(ON_pop_count, ON_case_count, by = c('Age', 'Gender'), all.x=F, all.y=F)
join_count = mutate(join_count, infection_rate = Freq/Total * 1000)



library(plotly)

Sys.setenv("plotly_username"="ags10001111")
Sys.setenv("plotly_api_key"="QvW27FIEknWIs4MLkZdu")

## Barplot of patients's condition in Ontario

p <- ggplot(join_count, aes(x = Age, y = infection_rate, fill = Gender)) +
  geom_col(position = 'dodge') + xlab('Age(years)') + ylab('Infection rate(%)') +
  ggtitle('Infection rate in Ontario(In percentage)') +
  geom_text(aes(label = paste(round(infection_rate, 2), '%')), position = position_dodge(0.9), vjust = -2, size = 4) +
  theme(text = element_text(size=10), 
        plot.title = element_text(hjust = 0.5, size = 20))

#p <- ggplotly(p)

levels(join_count$Gender) = c('Female', 'Male')

if(F){
p1 = plot_ly(join_count, x = ~Age, y = ~infection_rate, 
             type = "bar", color = ~ Gender) %>% 
             layout(title = 'Number of people infected per 1000 people in Ontario', 
                    xaxis = list(title = "Age (years)", tickangle = -45), 
                    yaxis = list(title = "Infected number"))
}

if(T){
  p1 = plot_ly(join_count, x = ~Age, y = ~infection_rate, 
               type = "bar", color = ~ Gender) %>% 
    layout(title = 'Infected cases (per 1,000 population) in Ontario', 
           xaxis = list(title = "Age (years)", tickangle = -45), 
           yaxis = list(title = "Infected cases"))
}

if(F){
  p1 = plot_ly(join_count, x = ~Age, y = ~infection_rate, 
               type = "bar", color = ~ Gender) %>% 
    layout(title = 'Infection rate in Ontario', 
           xaxis = list(title = "Age (years)", tickangle = -45), 
           yaxis = list(title = "Infection rate", tickformat = ".2%"))
}


# Choose a specific name for the plots.
api_create(p1, filename = "infection rate in ON")

