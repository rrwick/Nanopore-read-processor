#!/usr/bin/env Rscript

library(ggplot2)
library(scales)
library(tools)
library(grid)


normal_basecall_colour = "#2159aa" 
nanonet_basecall_colour = "#db7726"
very_good_to_bad_colours = c("#0571b0", "#92c5de", "#f28383", "#87001a")

# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      suppressWarnings(print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                              layout.pos.col = matchidx$col)))
    }
  }
}


# Get the input table from the command line argument.
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("You must supply one argument: the tsv file to analyse", call.=FALSE)
}
input_tsv = args[1]
read_data <- read.delim(input_tsv)
alignment_data_exists = "Alignment.identity" %in% colnames(read_data)


# Reorder the factors
read_data$Basecalling <- factor(read_data$Basecalling, levels = c("normal", "nanonet", "none"))
read_data$Quality.group <- factor(read_data$Quality.group, levels = c("very good", "good", "poor", "bad"))


# Make some subsets of the data
with_basecalling <- read_data[read_data$Basecalling != 'none' & !is.na(read_data$Basecalling),]
if (alignment_data_exists) {
  with_alignments <- with_basecalling[!is.na(with_basecalling$Alignment.identity),]
}


# Choose a good base count scale
total_base_count <- sum(with_basecalling$Length, na.rm = TRUE)
if (total_base_count > 2000000000) {
  label.format = scales::unit_format("Gb", 1e-9)
} else if (total_base_count > 2000000) {
  label.format = scales::unit_format("Mb", 1e-6)
} else {
  label.format = scales::unit_format("kb", 1e-3)
}

# Prepare some theme stuff.
my_theme <- theme_bw() + theme(aspect.ratio=1, plot.margin=margin(10, 10, 10, 10))
blank_theme <- my_theme + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank())
options(scipen=10000)


# Use smaller and more transparent points in the scatter plots when there are more reads.
point_size = 1 / (1 + (nrow(with_basecalling) / 2000))
point_alpha = 0.1 + (0.9 / (1 + (nrow(with_basecalling) / 10000)))






# Normal vs nanonet bases
normal_basecalls <- sum(read_data[read_data$Basecalling == "normal",]$Length, na.rm = TRUE)
nanonet_basecalls <- sum(read_data[read_data$Basecalling == "nanonet",]$Length, na.rm = TRUE)
bases <- data.frame(Basecalling = c("normal", "nanonet"),
                    value = c(normal_basecalls, nanonet_basecalls))
bases$Basecalling <- factor(bases$Basecalling, levels = c("normal", "nanonet"))
p1 <- ggplot(bases, aes(x=Basecalling, y=value, fill=Basecalling)) +
  ggtitle("Base count by basecalling method") + 
  geom_bar(colour = "black", stat = "identity") +
  my_theme +
  scale_y_continuous(labels=label.format) + 
  scale_fill_manual(values=c(normal_basecall_colour, nanonet_basecall_colour)) +
  guides(fill=FALSE) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())







# Length histogram
ninety_fifth_percentile <- quantile(with_basecalling$Length, 0.95)[[1]]
p2 <- ggplot(with_basecalling, aes(Length, fill = Basecalling)) +
  ggtitle("Length distribution") + 
  geom_histogram(binwidth = 50) +
  my_theme +
  scale_x_continuous(name="Read length (bp)") +
  scale_y_continuous(name="Read count") +
  coord_cartesian(xlim = c(0, ninety_fifth_percentile)) +
  scale_fill_manual(values=c(normal_basecall_colour, nanonet_basecall_colour)) +
  guides(fill=FALSE)




# Length histogram (log scale)
p3 <- ggplot(with_basecalling, aes(Length, fill = Basecalling)) +
  ggtitle("Length distribution (log scale)") + 
  geom_histogram(binwidth = 0.025) +
  my_theme +
  scale_x_continuous(name="Read length (bp)", trans = "log10", breaks = c(10, 100, 1000, 10000, 100000)) +
  scale_y_continuous(name="Read count") +
  coord_cartesian(xlim = c(10, 100000)) +
  scale_fill_manual(values=c(normal_basecall_colour, nanonet_basecall_colour)) +
  guides(fill=FALSE)






# These plots only apply if we have aligned our reads to a reference.
if (alignment_data_exists) {

  # Identity histogram
  p4 <- ggplot(with_alignments, aes(Alignment.identity, fill = Basecalling)) +
    ggtitle("Identity distribution") + 
    geom_histogram(binwidth = 0.5) +
    my_theme +
    scale_x_continuous(name="Alignment identity (%)", limits = c(45, 100), breaks = seq(40, 100, by = 5)) + 
    scale_y_continuous(name="Read count") +
    scale_fill_manual(values=c(normal_basecall_colour, nanonet_basecall_colour)) +
    guides(fill=FALSE)
  

  
  p5 <- ggplot(with_alignments, aes(x=Length, y=Alignment.identity, colour=Basecalling)) +
    ggtitle("Identity vs length") + 
    geom_point(alpha=point_alpha, size=point_size, shape=19) +
    my_theme +
    scale_x_continuous(name="Read length (bp)", trans = "log10", breaks = c(10, 100, 1000, 10000, 100000)) +
    scale_y_continuous(name="Alignment identity (%)") +
    coord_cartesian(xlim = c(10, 100000), ylim = c(45, 100)) +
    scale_colour_manual(values=c(normal_basecall_colour, nanonet_basecall_colour)) +
    guides(colour=FALSE)
  
  
  
  
  # Base quality quantity
  bad_bases <- sum(read_data[read_data$Quality.group == "bad",]$Length, na.rm = TRUE)
  poor_bases <- sum(read_data[read_data$Quality.group == "poor",]$Length, na.rm = TRUE)
  good_bases <- sum(read_data[read_data$Quality.group == "good",]$Length, na.rm = TRUE)
  very_good_bases <- sum(read_data[read_data$Quality.group == "very good",]$Length, na.rm = TRUE)
  bases <- data.frame(Quality = c("very good", "good", "poor", "bad"),
                      value = c(very_good_bases, good_bases, poor_bases, bad_bases))
  bases$Quality <- factor(bases$Quality, levels = c("very good", "good", "poor", "bad"))
  p6 <- ggplot(bases, aes(x=Quality, y=value, fill=Quality)) +
    ggtitle("Base count by quality") + 
    geom_bar(colour = "black", stat = "identity") +
    my_theme +
    scale_y_continuous(labels=label.format) +
    scale_fill_manual(values=very_good_to_bad_colours) +
    guides(fill=FALSE) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
  
  
  p7 <- ggplot(with_alignments, aes(x=Length, y=Alignment.identity, colour=Quality.group)) +
    ggtitle("Identity vs length") + 
    geom_point(alpha=point_alpha, size=point_size, shape=19) +
    my_theme +
    scale_x_continuous(name="Read length (bp)", trans = "log10", breaks = c(10, 100, 1000, 10000, 100000)) +
    scale_y_continuous(name="Alignment identity (%)") +
    coord_cartesian(xlim = c(10, 100000), ylim = c(45, 100)) +
    scale_colour_manual(values=very_good_to_bad_colours) +
    guides(colour=FALSE)
  
  
  
  # Qscore vs identity scatter
  p8 <- ggplot(with_alignments, aes(x=Mean.qscore, y=Alignment.identity, colour=Basecalling)) +
    ggtitle("Identity vs qscore") + 
    geom_point(alpha=point_alpha, size=point_size, shape=19) +
    my_theme +
    scale_x_continuous(name="Mean qscore") +
    scale_y_continuous(name="Alignment identity (%)") +
    coord_cartesian(xlim = c(0, 20), ylim = c(45, 100)) +
    scale_colour_manual(values=c(normal_basecall_colour, nanonet_basecall_colour)) +
    guides(colour=FALSE)
}


# Save plots to png
png_name = paste(file_path_sans_ext(input_tsv), "_plots.png", sep="")
options(bitmapType='cairo')
if (alignment_data_exists) {
  png(png_name, width = 4000, height = 4000, res = 300)
  multiplot(p1, p4, p6, p2, p5, p7, p3, p8, cols=3)
} else {
  png(png_name, width = 4000, height = 1200, res = 300)
  multiplot(p1, p2, p3, cols=3)
}
garbage <- dev.off()


# Find the n50 read length
sorted_read_lengths = sort(with_basecalling$Length, decreasing = TRUE)
n50_target = sum(with_basecalling$Length) / 2
x <- 1
while(sum(sorted_read_lengths[1:x]) < n50_target) {x <- x+1;}
n50 = sorted_read_lengths[x]


# Print some statistics
message(paste("Mean read length:   ", round(mean(with_basecalling$Length), 2), sep=""))
message(paste("Median read length: ", median(with_basecalling$Length), sep=""))
message(paste("N50 read length:    ", n50, sep=""))
if (alignment_data_exists) {
  message(paste("Mean identity:      ", round(mean(with_alignments$Alignment.identity), 2), "%", sep=""))
  message(paste("Median identity:    ", median(with_alignments$Alignment.identity), "%", sep=""))
}
