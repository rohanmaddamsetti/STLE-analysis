## windows.R by Rohan Maddamsetti.

library(ggplot2)

## Make a histogram of the size of resident regions that survived the recombination
## treatment in the Souza-Turner experiment.

window.data <- read.csv("big_windows_REL4397.csv")
window.data$Length <- window.data$End - window.data$Start

quartz()
window.hist <- ggplot(window.data,aes(x=Length)) + geom_histogram() #+ scale_x_log10()
window.hist
