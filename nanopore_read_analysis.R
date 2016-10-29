#!/usr/bin/env Rscript

library(ggplot2)
library(scales)
library(tools)
library(grid)

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


# Reorder the Basecalling factors
read_data$Basecalling <- factor(read_data$Basecalling, levels = c("normal", "nanonet", "none"))


# Make some subsets of the data
with_basecalling <- read_data[read_data$Basecalling != 'none',]
normal_basecalling = read_data[read_data$Basecalling == "normal",]
nanonet_basecalling = read_data[read_data$Basecalling == "nanonet",]
if (alignment_data_exists) {
  with_alignments <- with_basecalling[!is.na(with_basecalling$Alignment.identity),]
}

# Prepare some theme stuff.
my_theme <- theme_bw() + theme(aspect.ratio=1, plot.margin=margin(10, 10, 10, 10))
blank_theme <- my_theme + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank())
options(scipen=10000)


# Use smaller and more transparent points in the scatter plots when there are more reads.
point_size = 1 / (1 + (nrow(with_basecalling) / 2000))
point_alpha = 0.1 + (0.9 / (1 + (nrow(with_basecalling) / 10000)))


# Length histogram
ninety_ninth_percentile <- quantile(with_basecalling$Length, 0.99)[[1]]
p1 <- ggplot(with_basecalling, aes(Length, fill = Basecalling)) +
  ggtitle("Length distribution") + 
  geom_histogram(binwidth = 100) +
  my_theme +
  scale_x_continuous(name="Read length (bp)", limits = c(0, ninety_ninth_percentile)) +
  scale_y_continuous(name="Read count")


# Length histogram (log scale)
p2 <- ggplot(with_basecalling, aes(Length, fill = Basecalling)) +
  ggtitle("Length distribution (log scale)") + 
  geom_histogram(binwidth = 0.025) +
  my_theme +
  scale_x_continuous(name="Read length (bp)", trans = "log10", limits = c(100, 100000), breaks = c(100, 1000, 10000, 100000)) +
  scale_y_continuous(name="Read count")


# Total bases pie chart
normal_basecalls <- sum(normal_basecalling$Length)
nanonet_basecalls <- sum(nanonet_basecalling$Length)
bases <- data.frame(Basecalling = c("normal", "nanonet"),
                    value = c(normal_basecalls, nanonet_basecalls),
                    mid_y = c(nanonet_basecalls / 2, nanonet_basecalls + (normal_basecalls / 2)))
bases$Basecalling <- factor(bases$Basecalling, levels = c("normal", "nanonet"))
text = c(paste(prettyNum(nanonet_basecalls, big.mark=",", scientific=FALSE, preserve.width="none"), " bp"),
         paste(prettyNum(normal_basecalls, big.mark=",", scientific=FALSE, preserve.width="none"), " bp"))
p3 <- ggplot(bases, aes(x="", y=value, fill=Basecalling)) +
  ggtitle("Total bases") + 
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + 
  blank_theme +
  theme(axis.text.x=element_blank()) + 
  geom_text(aes(y = mid_y, label = text))


# These plots only apply if we have aligned our reads to a reference.
if (alignment_data_exists) {

  # Identity histogram
  p4 <- ggplot(with_alignments, aes(Alignment.identity, fill = Basecalling)) +
    ggtitle("Identity distribution") + 
    geom_histogram(binwidth = 0.5) +
    my_theme +
    scale_x_continuous(name="Alignment identity (%)", limits = c(45, 100), breaks = seq(40, 100, by = 5)) + 
    scale_y_continuous(name="Read count")
  
  
  # Length vs identity scatter
  p5 <- ggplot(with_alignments, aes(x=Length, y=Alignment.identity, colour=Basecalling)) +
    ggtitle("Identity against length") + 
    geom_point(alpha=point_alpha, size=point_size, shape=19) +
    my_theme +
    guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
    scale_x_continuous(name="Read length (bp)", trans = "log10", limits = c(100, 100000), breaks = c(100, 1000, 10000, 100000)) +
    scale_y_continuous(name="Alignment identity (%)", limits = c(50, 100))
  
  
  # Qscore vs identity scatter
  p6 <- ggplot(with_alignments, aes(x=Mean.qscore, y=Alignment.identity, colour=Basecalling)) +
    ggtitle("Identity against qscore") + 
    geom_point(alpha=point_alpha, size=point_size, shape=19) +
    my_theme +
    guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
    scale_x_continuous(name="Mean qscore", limits = c(0, 20)) +
    scale_y_continuous(name="Alignment identity (%)", limits = c(50, 100))
}


# Save plots to png
png_name = paste(file_path_sans_ext(input_tsv), "_plots.png", sep="")
options(bitmapType='cairo')
if (alignment_data_exists) {
  png(png_name, width = 5000, height = 3000, res = 300)
  multiplot(p1, p4, p2, p5, p3, p6, cols=3)
} else {
  png(png_name, width = 5000, height = 1500, res = 300)
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
message(paste("Mean identity:      ", round(mean(with_alignments$Alignment.identity), 2), "%", sep=""))
message(paste("Median identity:    ", median(with_alignments$Alignment.identity), "%", sep=""))
