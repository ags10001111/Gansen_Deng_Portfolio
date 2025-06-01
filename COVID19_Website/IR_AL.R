setwd("C:/Canteen/Western/2020 Summer/COVID19/Plot/Infection rate")
library(tidyverse)

## Import and attach the data for Ontario
mydat <- read.csv('covid19dataexport.csv', header=TRUE)
mydat = mydat[mydat$Gender == 'Male' | mydat$Gender == 'Female', ]
AL_case = mydat[, which(colnames(mydat) %in% c('Gender', 'Age.group'))]

## Remove unused levels
AL_case = AL_case[AL_case$Age.group != 'Unknown', ]
AL_case = droplevels(AL_case)

## Read Alberta population data
AL_pop = read.csv('AL_pop.csv')
AL_pop$Age = AL_pop$ï..Age
AL_pop = droplevels(AL_pop)


## Reassign levels
levels(AL_case$Age.group) = c(rep('0 - 19 ', 2), '20 - 29 ', '30 - 39 ',
                              '40 - 49 ', '0 - 19 ', '50 - 59 ', '60 - 69 ', 
                              '70 - 79 ', 'over 80 ', '0 - 19 ')


levels(AL_pop$Age) = c(rep('0 - 19 ', 2), 'over 80 ', '0 - 19 ',
                       rep('20 - 29 ', 2), rep('30 - 39 ', 2), 
                       rep('40 - 49 ', 2), '0 - 19 ', 
                       rep('50 - 59 ', 2), 
                       rep('60 - 69 ', 2), rep('70 - 79 ', 2), 
                       rep('over 80 ', 4))

# Remove commas from the number
library(tidyverse)
AL_pop$Count = as.numeric(gsub(",","",AL_pop$Count))
AL_pop_count = AL_pop %>% 
  group_by(Age, Gender) %>% 
  summarise(Total = sum(Count))

AL_case_count = as.data.frame(table(AL_case))
colnames(AL_case_count) = c('Gender', 'Age', 'Freq')
join_count = merge(AL_pop_count, AL_case_count, by = c('Age', 'Gender'), all.x=F, all.y=F)
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

p <- ggplotly(p)

join_count$Age = factor(join_count$Age, levels(join_count$Age)[c(1, 3:8, 2)])

if(T){
p1 = plot_ly(join_count, x = ~Age, y = ~infection_rate, 
             type = "bar", color = ~ Gender) %>% 
  layout(title = 'Infected cases (per 1,000 population) in Alberta', 
         xaxis = list(title = "Age (years)", tickangle = -45), 
         yaxis = list(title = "Infected cases", tickformat = ".2f")) 
}

if(F){
  p1 = plot_ly(join_count, x = ~Age, y = ~infection_rate, 
               type = "bar", color = ~ Gender) %>% 
    layout(title = 'Infection rate in Alberta', 
           xaxis = list(title = "Age (years)", tickangle = -45), 
           yaxis = list(title = "Infection rate", tickformat = ".2%")) 
}


# Choose a specific name for the plots.
api_create(p1, filename = "infection rate in AL")
