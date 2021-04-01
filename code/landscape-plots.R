library("tidyverse")
library("ggtree")
library("RColorBrewer")

# wipe memory and construct list of files
rm(list=ls())
files <- commandArgs(trailingOnly = T)
prefix <- "results/landscape-tables-without-CpGmod"
figures.dir <- 'results/figures'
treefile <- "doc/taxon-sampling.tree"
landscapes.dir <- file.path(figures.dir, 'landscapes-without-CpGmod')

if (length(files) == 0) {
	files <- list.files(path = prefix, pattern = "*.txt")
	files <- file.path(prefix, files)
}



# prepare data frame
divcol <- c()
typecol <- c()
valuecol <- c()
speciescol <- c()
tablata = data.frame(speciescol, divcol, typecol, valuecol)

# load all the tables to construct a single data frame
for (thisfile in files) {
	# species name needs to be in the file name
	species <- sub(".txt$", "", basename(thisfile))
	species <- sub("_", " ", species)
	writeLines(paste("Loading", thisfile, "for", species))
	data <- read.table(thisfile, header = T)

	# reformat the data frame so it looks like this:
	# Divergence    Type     Value
	# 0          Unknown 0.6970246
	# 1          Unknown 0.9783817

	# get the list of types
	types <- names(data)[2:length(names(data))]

	# get the list of divergence levels
	divergence <- data[,1]

	# generate the three columns:
	# divergence must simply be replicated
	divcol <- rep(divergence, length(types))
	# the other two are filled from the data frame
	typecol <- c()
	valuecol <- c()
	for (type in types) {
		# the type must simply be replicated 51 times
		typecol <- c(typecol, rep(type, 51))
		# the values are in the table in the type column
		valuecol <- c(valuecol, data[,type])
	}
	speciescol <- rep(species, length(types))

	# construct new data frame
	temptable <- data.frame(speciescol, divcol, typecol, valuecol)
	names(temptable) <- c("Species", "Divergence", "Type", "Value")
	tablata <- rbind(tablata, temptable)
}


# species sorted by phylogenetic placement, required for species order
tree <- read.tree(treefile)
tree$tip.label <- sub("_", " ", tree$tip.label) # replace _ with space
mybreaks <- rev(tree$tip.label)

# reorder the species levels
tablata <- within(tablata, Species <- factor(Species, levels = mybreaks))

# make the Divergence a real Kimura distance (between 0 and 1)
tablata$Divergence <- tablata$Divergence / 100

writeLines(paste('Creating plots, they will all be saved in', figures.dir))
writeLines(c(rep('-', 44), '\n'), sep='')


getPalette <- function(x) {
	# first get a color palette that matches the number of stuff
	palette <- colorRampPalette(brewer.pal(9, "Set1"))
	colorcount <- length(levels(x))
	# and make it a named vector
	pal <- palette(colorcount)
	names(pal) <- sort(levels(x))
	return(pal)
}

# plot dat, make barcharts
label_y <- "Genome coverage [%]"
label_y_log <- "Genome coverage [log %]"
label_x <- "Kimura distance from TE family consensus sequence"
mypalette <- getPalette(tablata$Type)


# create plots for each species
# this takes long
for (thisspecies in levels(tablata$Species)) {
	writeLines(paste("Creating plots for", thisspecies))
	thisdata <- subset(tablata, Species == thisspecies)
	thistitle <- paste("Repeat landscape for", thisspecies)

	# bar chart
	barchart <- ggplot(thisdata, aes(x = Divergence, y = Value)) + geom_bar(stat="identity", aes(fill = Type))
	barchart <- barchart + scale_fill_manual(values = mypalette, breaks = names(mypalette)) + theme_bw()
	barchart <- barchart + xlab(label_x) + ylab(label_y) + ggtitle(thistitle)
	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-barchart.pdf")), plot = barchart, units="mm", width=297, height=210)
	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-barchart.svg")), plot = barchart, units="mm", width=297, height=210)

 	# heatmap
 	heatmap <- ggplot(thisdata, aes(x = Divergence, y = Type)) + geom_tile(aes(fill = Value))
 	heatmap <- heatmap + ylab("TE family") + xlab(label_x) + labs(fill=label_y)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-heatmap.pdf")), plot=heatmap, units="mm", width=210, height=297)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-heatmap.svg")), plot=heatmap, units="mm", width=210, height=297)

 	# heatmap again, with log scale
 	heatmap <- ggplot(thisdata, aes(x = Divergence, y = Type)) + geom_tile(aes(fill = log(Value)))
 	heatmap <- heatmap + ylab("TE family") + xlab(label_x) + labs(fill=label_y_log)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-heatmap-logscale.pdf")), plot=heatmap, units="mm", width=210, height=297)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-heatmap-logscale.svg")), plot=heatmap, units="mm", width=210, height=297)

 	# bar chart again, without unknown
 	writeLines("Removing 'Unknown' elements from the data set and creating another bar chart")
 	thisdata <- subset(thisdata, Type != 'Unknown')
 	barchart <- ggplot(thisdata, aes(x = Divergence, y = Value)) + geom_bar(stat="identity", aes(fill = Type))
 	barchart <- barchart + scale_fill_manual(values = mypalette, breaks = names(mypalette)) + theme_bw()
 	barchart <- barchart + xlab(label_x) + ylab(label_y) + ggtitle(thistitle)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-nounknown-barchart.pdf")), plot = barchart, units="mm", width=297, height=210)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-nounknown-barchart.svg")), plot = barchart, units="mm", width=297, height=210)
 	ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-nounknown-barchart.png")), plot = barchart, units="mm", width=210, height=148)

}

# reduce the TE types to class/superfamilies DNA, LTR, LINE, SINE
tablata$Type <- sub("\\..+", "", tablata$Type)
tablata$Type <- sub("Retroposon", "LTR", tablata$Type)
mylimits <- c("DNA", "LTR", "LINE", "SINE", "RC", "Other", "Unknown")
tablata$Type <- ordered(tablata$Type, levels = mylimits)
# new palette for reduced levels
mypalette <- getPalette(tablata$Type)
names(mypalette) <- mylimits
label_y <- "Genome coverage [%]"
label_x <- "Kimura distance from TE family consensus sequence"

# this saves a faceted barplot with all species (awesome!) :)
writeLines("Creating bar chart with all species")
title <- ""

# construct the bar chart first -- this by itself cannot be printed as the displayed data are summarized
barchart <- ggplot(tablata, aes(x = Divergence, y = Value))
barchart <- barchart + geom_bar(stat="identity", aes(fill = Type))
barchart <- barchart + scale_fill_manual(values = mypalette, breaks = names(mypalette))
barchart <- barchart + xlab(label_x) + ylab(label_y) + ggtitle(title)
barchart <- barchart + theme_bw() + theme(strip.text = element_text(face = "italic"))

# make the posters by using facet_wrap by species
poster <- barchart + facet_wrap(~ Species, scales = "free_y") + guides(fill = guide_legend(ncol=1)) # legend in two columns because too long
ggsave(file.path(figures.dir, "landscape-all-A3.pdf"), plot=poster, units="mm", width=420, height=297) # DIN A3 landscape
ggsave(file.path(figures.dir, "landscape-all-A3.svg"), plot=poster, units="mm", width=420, height=297) # DIN A3 landscape, svg
poster <- barchart + facet_wrap(~ Species, scales = "free_y") + guides(fill = guide_legend(ncol=1)) # legend in single column for portrait mode
ggsave(file.path(figures.dir, "landscape-all-A2.pdf"), plot=poster, units="mm", width=594, height=420) # DIN A2 landscape
ggsave(file.path(figures.dir, "landscape-all-A2.svg"), plot=poster, units="mm", width=594, height=420) # DIN A2 landscape, svg
ggsave(file.path(figures.dir, "landscape-all-A1.pdf"), plot=poster, units="mm", width=594, height=841) # DIN A1 portrait
ggsave(file.path(figures.dir, "landscape-all-A1.svg"), plot=poster, units="mm", width=594, height=841) # DIN A1 portrait, svg

# remove unknown elements and do the same again
writeLines("Creating bar chart with all species and without Unknown elements")
tablata.nounknown <- subset(tablata, Type != 'Unknown')

barchart <- ggplot(tablata.nounknown, aes(x = Divergence, y = Value)) + geom_bar(stat="identity", aes(fill = Type))
barchart <- barchart + scale_fill_manual(values = mypalette, breaks = names(mypalette)) + xlab(label_x) + ylab(label_y) + ggtitle(title) + theme_bw() + theme(strip.text = element_text(face = "italic"))
poster <- barchart + facet_wrap(~ Species, scales = "free_y") + guides(fill = guide_legend(ncol=1)) # legend in two columns because too long
ggsave(file.path(figures.dir, "landscape-all-nounknown-A3.pdf"), plot=poster, units="mm", width=420, height=297) # DIN A3 landscape
ggsave(file.path(figures.dir, "landscape-all-nounknown-A3.svg"), plot=poster, units="mm", width=420, height=297) # DIN A3 landscape, svg
poster <- barchart + facet_wrap(~ Species, scales = "free_y") + guides(fill = guide_legend(ncol=1)) # legend in single column for portrait mode
ggsave(file.path(figures.dir, "landscape-all-nounknown-A2.pdf"), plot=poster, units="mm", width=594, height=420) # DIN A2 landscape
ggsave(file.path(figures.dir, "landscape-all-nounknown-A2.svg"), plot=poster, units="mm", width=594, height=420) # DIN A2 landscape, svg
ggsave(file.path(figures.dir, "landscape-all-nounknown-A1.pdf"), plot=poster, units="mm", width=594, height=841) # DIN A1 portrait
ggsave(file.path(figures.dir, "landscape-all-nounknown-A1.svg"), plot=poster, units="mm", width=594, height=841) # DIN A1 portrait, svg format



# I want to add the landscape plots to the tree tips, so let's make a list of ggplot objects
bx <- list()
for (thisspecies in tree$tip.label) {
	writeLines(paste("Creating reduced plots for", thisspecies))
	thisdata <- subset(tablata.nounknown, Species == thisspecies)

	# bar chart
	barchart <- ggplot(thisdata, aes(x = Divergence, y = Value)) + geom_bar(stat="identity", aes(fill = Type))
	barchart <- barchart + scale_fill_manual(values = mypalette, breaks = names(mypalette)) + theme_inset()
	barchart <- barchart + theme(axis.title = element_blank())
	barchart <- barchart + theme(legend.position = "none")
	barchart <- barchart + theme(axis.text = element_text(size = 16))
	barchart <- barchart + theme(axis.ticks = element_line(size = 1))
	bx[[thisspecies]] <- barchart
	# DIN A5 (ISO 216)
	# ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-reduced-nounknown-barchart.pdf")), plot = barchart, units="mm", width=210, height=148)
	# ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-reduced-nounknown-barchart.svg")), plot = barchart, units="mm", width=210, height=148)
	# ggsave(file.path(landscapes.dir, paste0(thisspecies, "-repeat-landscape-reduced-nounknown-barchart.png")), plot = barchart, units="mm", width=210, height=80)
}
names(bx) <- 1:73 # need to number the list elements because ggtree uses node IDs for tips, not their labels

# plot settings
fontsize <- 5
offset <- 9
offset.text <- 0.3

# make tree plot
p <- ggtree(tree) + geom_tiplab(size = fontsize) + xlim(0,40)
# add clade labels
p <- p +
	geom_cladelabel(node=87,    offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#343d46", fontsize = fontsize, label = "Diptera") +
	geom_cladelabel(node = 107, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#4f5b66", fontsize = fontsize, label = "Lepidoptera") +
	geom_cladelabel(node = 112, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#65737e", fontsize = fontsize, label = "Coleoptera") +
	geom_cladelabel(node = 116, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#536872", fontsize = fontsize, label = "Hymenoptera") +
	geom_cladelabel(node = 129, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#343d46", fontsize = fontsize, label = "Hemiptera") +
	geom_cladelabel(node = 138, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#4f5b66", fontsize = fontsize, label = "Palaeoptera") +
	geom_cladelabel(node = 139, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#65737e", fontsize = fontsize, label = "Crustacea") +
	geom_cladelabel(node = 141, offset = offset, angle = 90, offset.text = offset.text, hjust = 0.5, barsize = 2, color = "#536872", fontsize = fontsize, label = "Chelicerata")
# add landscape plots
treeplot <- inset(p, bx, height = 0.017, width = 0.3, hjust = -14)
# tadaa
ggsave("results/figures/tree-with-landscapes.pdf", plot = treeplot, units = "mm", width = 294, height = 841)
ggsave("results/figures/tree-with-landscapes.svg", plot = treeplot, units = "mm", width = 294, height = 841)
ggsave("results/figures/tree-with-landscapes.png", plot = treeplot, units = "mm", width = 294, height = 841, dpi = 300)
