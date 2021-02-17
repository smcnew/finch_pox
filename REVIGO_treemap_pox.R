

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0009607","response to biotic stimulus",0.342,2.3625,0.787,0.000,"response to biotic stimulus"),
c("GO:0009605","response to external stimulus",1.370,2.0448,0.770,0.370,"response to biotic stimulus"),
c("GO:0035634","response to stilbenoid",0.001,2.2733,0.829,0.218,"response to biotic stimulus"),
c("GO:0032480","negative regulation of type I interferon production",0.004,3.9469,0.500,0.000,"negative regulation of type I interferon production"),
c("GO:1902957","negative regulation of mitochondrial electron transport, NADH to ubiquinone",0.000,2.6696,0.566,0.332,"negative regulation of type I interferon production"),
c("GO:1902524","positive regulation of protein K48-linked ubiquitination",0.000,2.6696,0.753,0.118,"negative regulation of type I interferon production"),
c("GO:0045824","negative regulation of innate immune response",0.012,3.4609,0.453,0.399,"negative regulation of type I interferon production"),
c("GO:0044033","multi-organism metabolic process",0.025,2.3696,0.876,0.000,"multi-organism metabolism"),
c("GO:0000270","peptidoglycan metabolic process",0.769,2.0701,0.808,0.492,"multi-organism metabolism"),
c("GO:0051672","catabolism by organism of cell wall peptidoglycan in other organism",0.000,2.4949,0.508,0.425,"multi-organism metabolism"),
c("GO:0071554","cell wall organization or biogenesis",0.950,2.1278,0.945,0.000,"cell wall organization or biogenesis"),
c("GO:0019835","cytolysis",0.044,2.1278,0.949,0.026,"cytolysis"),
c("GO:0034354","'de novo' NAD biosynthetic process from tryptophan",0.026,2.1945,0.709,0.044,"primede novoprime NAD biosynthesis from tryptophan"),
c("GO:0044036","cell wall macromolecule metabolic process",0.709,2.1278,0.897,0.045,"cell wall macromolecule metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_pox.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
