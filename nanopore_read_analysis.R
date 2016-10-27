library(ggplot2)
library(scales)

read_data <- read.delim("~/Desktop/Kleb_INF167_2d.tsv")

# Reorder the Basecalling factors
read_data$Basecalling <- factor(read_data$Basecalling, levels = c("normal", "nanonet", "none"))

# Make some subsets of the data
with_basecalling <- read_data[read_data$Basecalling != 'none',]
normal_basecalling = read_data[read_data$Basecalling == "normal",]
nanonet_basecalling = read_data[read_data$Basecalling == "nanonet",]

blank_theme <- theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
options(scipen=10000)

# Length histogram
ggplot(with_basecalling, aes(Length, fill = Basecalling)) +
  ggtitle("Read length distribution") + 
  geom_histogram(binwidth = 100) +
  theme_bw() +
  scale_x_continuous(name="Read length (bp)", limits = c(0, 20000)) +
  scale_y_continuous(name="Read count")

# Length histogram (log scale)
ggplot(with_basecalling, aes(Length, fill = Basecalling)) +
  ggtitle("Read length distribution (log)") + 
  geom_histogram(binwidth = 0.025) +
  theme_bw() +
  scale_x_continuous(name="Read length (bp)", trans = "log10", limits = c(100, 100000), breaks = c(100, 1000, 10000, 100000)) +
  scale_y_continuous(name="Read count")

# Identity histogram
ggplot(with_basecalling, aes(Alignment.identity, fill = Basecalling)) +
  ggtitle("Alignment identity distribution") + 
  geom_histogram(binwidth = 0.5) +
  theme_bw() +
  scale_x_continuous(name="Alignment identity (%)", limits = c(45, 100), breaks = seq(40, 100, by = 5)) + 
  scale_y_continuous(name="Read count")

# Length vs identity scatter
ggplot(with_basecalling, aes(x=Length, y=Alignment.identity, colour=Basecalling)) +
  ggtitle("Read length vs alignment identity") + 
  geom_point(alpha=0.2, size=0.25) +
  theme_bw() +
  scale_x_continuous(name="Read length (bp)", trans = "log10", limits = c(100, 100000), breaks = c(100, 1000, 10000, 100000)) +
  scale_y_continuous(name="Alignment identity (%)", limits = c(50, 100))

# Qscore vs identity scatter
ggplot(with_basecalling, aes(x=Mean.qscore, y=Alignment.identity, colour=Basecalling)) +
  ggtitle("Mean qscore vs alignment identity") + 
  geom_point(alpha=0.2, size=0.25) +
  theme_bw() +
  scale_x_continuous(name="Mean qscore", limits = c(0, 20)) +
  scale_y_continuous(name="Alignment identity (%)", limits = c(50, 100))

# Total bases - normal and nanonet basecalling
bases <- data.frame(Basecalling = c("normal", "nanonet"),
                    value = c(sum(normal_basecalling$Length), sum(nanonet_basecalling$Length)))
ggplot(bases, aes(x="", y=value, fill=Basecalling)) +
  ggtitle("Total bases") + 
  geom_bar(width = 1, stat = "identity") + 
  theme_bw() + blank_theme + 
  theme(axis.text.x=element_blank()) + 
  coord_polar("y", start=0) +
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                label = paste(prettyNum(value, big.mark=",", scientific=FALSE, preserve.width="none"), " bp")),
            size=5)
