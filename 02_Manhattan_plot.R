########################Manhattan_plot########################
#This R script was used to generate and save Manhattan plots for 
#psoriasis and PRS phewas, stratified by sex.

#by: Luo wenjin
#date: 2024-01-15
##############################################################


library(data.table) 
library(ggplot2)     
library(ggrepel)     
library(dplyr)       

# Set working directory
setwd('setwd("~/data/PSO/manhattan_plot/")  

# Define chromosome labels
chr_lab <- c('lip, oral cavity and pharynx',
             'digestive organs',
             'respiratory organs',
             'bone and articular cartilage',
             'skin',
             'mesothelial and soft tissue',
             'breast',
             'female genital organs',
             'male genital organs',
             'urinary tract',
             'central nervous system',
             'endocrine glands',
             'ISUS',
             'haematopoietic',
             'primary multiple sites')

# Define function to create Manhattan plot
manhattan_plot <- function(x, y) {
  dt1 <- x  # Input data
  pos <- dt1$pvalue  # Extract p-values
  counters <- rep(0, length(pos))  # Initialize counters
  count <- 0
  
  # Count consecutive non-NA p-values
  for (i in 1:length(pos)) {
    if (is.na(pos[i])) {
      count <- 0
    } else {
      count <- count + 1
      counters[i] <- count
    }
  }
  
  num <- c(counters[which(is.na(pos)) - 1], 1)  # Count non-NA groups
  chr <- vector()  # Initialize chromosome vector
  
  # Assign chromosome numbers based on group sizes
  for (i in 1:length(num)) {
    chr <- append(chr, rep(i, num[i]))
  }
  
  # Prepare data for plotting
  dt2 <- dt1 %>%
    filter(!is.na(pvalue)) %>%
    mutate('CHR' = chr)
  
  n_row <- nrow(dt2)  # Number of rows in filtered data
  dt2 <- dt2 %>% mutate('bp' = seq(1, n_row * 3, 3))  # Assign base pair positions
  
  tb1 <- dt2
  
  # Calculate cumulative base pair positions
  data_cum <- tb1 %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(max_bp = max(bp)) %>%
    dplyr::mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(CHR, bp_add)
  
  tb1 <- tb1 %>%
    dplyr::inner_join(data_cum, by = "CHR") %>%
    dplyr::mutate(bp_cum = bp + bp_add)
  
  # Set axis labels
  axis_set <- tb1 %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarize(center = mean(bp_cum))
  
  # Define y-axis limits
  ylim <- tb1 %>%
    filter(pvalue == min(pvalue)) %>%
    mutate(ylim = abs(floor(log10(pvalue))) + 2) %>%
    pull(ylim)
  
  # Create Manhattan plot
  p1 <- ggplot(tb1, aes(x = bp_cum, y = -log10(pvalue), color = as.factor(CHR))) +
    geom_hline(yintercept = seq(0, 4, 1), color = "#DCDCDC", 
               linetype = "solid", alpha = 1) + 
    geom_point(data = tb1, alpha = 1, shape = 16, size = 7) +
    geom_text_repel(data = tb1[tb1$pvalue < 0.05,], 
                    aes(label = Cancer), 
                    color = "black") +
    scale_x_continuous(label = y, breaks = axis_set$center) +
    scale_y_continuous(breaks = seq(0, 6, 1)) + 
    scale_color_manual(values = rep(c("#91CCC0","#7FABD1","#EC6E66"), 
                                    unique(length(axis_set$CHR)))) +
    labs(x = NULL, y = "-log10(P)") +
    theme_test(base_size = 15,
               base_family = "arial") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1),
          legend.position = "none") 

  return(p1)  # Return the plot
}

# Read data from CSV files
df1 <- fread('plot1.csv')
df2 <- fread('plot2.csv')
df3 <- fread('plot3.csv')
df4 <- fread('plot4.csv')
df5 <- fread('plot5.csv')
df6 <- fread('plot6.csv')

# Generate Manhattan plots and save to PDF files
pdf("manhattan_plot_df1.pdf")
print(manhattan_plot(df1, chr_lab[-c(8,9)]))
dev.off()

pdf("manhattan_plot_df2.pdf")
print(manhattan_plot(df2, chr_lab[-8]))
dev.off()

pdf("manhattan_plot_df3.pdf")
print(manhattan_plot(df3, chr_lab[-9]))
dev.off()

pdf("manhattan_plot_df4.pdf")
print(manhattan_plot(df4, chr_lab[-c(8,9)]))
dev.off()

pdf("manhattan_plot_df5.pdf")
print(manhattan_plot(df5, chr_lab[-8]))
dev.off()

pdf("manhattan_plot_df6.pdf")
print(manhattan_plot(df6, chr_lab[-9]))
dev.off()
