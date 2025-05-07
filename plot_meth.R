dir = '/Users/dashnowh/Documents/git/inSTRbility'
setwd(dir)
library(ggplot2)
library(cowplot)
# apply cowplot theme
theme_set(theme_cowplot())

# Load the data into a single data frame
all.data = data.frame()
for (file in list.files(path = dir, pattern = '*meth.txt', full.names = TRUE)) {
  this.data = read.table(file, header = F, stringsAsFactors = F)
  this.sample = gsub('.meth.txt', '', basename(file))
  this.data$sample = this.sample
  all.data = rbind(all.data, this.data)

}

colnames(all.data) = c('repeatlen', 'medianmeth', 'sample')
all.data$medianmeth = as.numeric(all.data$medianmeth)
all.data$medianmethlevel = round(all.data$medianmeth, 1)
all.data$ismeth = all.data$medianmeth >= 0.5

ggplot(all.data, aes(x = repeatlen, fill = as.factor(medianmethlevel))) + 
  geom_histogram() + facet_wrap(~sample, scales = 'free') + 
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation')  
