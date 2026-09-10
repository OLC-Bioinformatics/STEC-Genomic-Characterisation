#if (!requireNamespace("BiocManager", quietly=TRUE))
# install.packages("BiocManager")
#BiocManager::install("BiocUpgrade") ## you may need this
#BiocManager::install("ggtreeExtra")
#BiocManager::install("YuLab-SMU/ggtree")

#install chatgpt
# install.packages("chattr")
# library(chattr)
#install.packages("janitor")


#now load the packages
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(TDbook)
#import the treefile
library(ape)
#rename the tree
library(treeio)
library(dplyr) #dataframe manipulation
library(ggnewscale)#to add a new color scale and not mess with original color scale
library(readxl) #to read xlsx files
library(paletteer) #colour palettes
paletteer_d("colorBlindness::Blue2Orange10Steps")
library(RColorBrewer) #color palettes
library(janitor) # package to clean up column names in table (make all lower case, underscores instead of spaces)

#the colorblind friendly brewer palettes:
display.brewer.all(colorblindFriendly = TRUE)

#set our working directory
setwd("/Users/noorshubair/Desktop/chapter4/chewbbaca/cgmlst_0.95/")

#import the treefile
nwk <- read.tree(file="cgmlst95_MSTreeV2.nwk")
nwk
class(nwk)

length(nwk$tip.label)
ggtree(nwk)


#lets write the tip labels to a csv file so I can create a metadatatable in excel
nwk$tip.label
write.csv(nwk$tip.label, file="nwktreetiplabels.csv")

#have a look at the tree (note the branch lengths are missing here, so it is round)
ggtree(nwk, 
       layout = "daylight", branch.length = "none") +
  theme_tree()

#bring in some metadata
nwkmeta <- read_excel(
  "/Users/noorshubair/Desktop/chapter4/data/Stx2i_genomes_working_260730_forR_cgmlstfigures.xlsx",
  sheet = "NCBI_CFIA_dataset"
)
nwkmeta

#clean up column names
nwkmeta <- nwkmeta %>%
  clean_names() %>%
  mutate(label = assembly) %>%
  relocate(label, .before = 1)

# Confirm that every tree tip has a matching metadata row
missing_metadata <- setdiff(nwk$tip.label, nwkmeta$label)
extra_metadata <- setdiff(nwkmeta$label, nwk$tip.label)

if (length(missing_metadata) > 0) {
  stop(
    "These tree tips are missing from the metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

if (length(extra_metadata) > 0) {
  warning(
    "These metadata rows are not present in the tree: ",
    paste(extra_metadata, collapse = ", ")
  )
}

g <- ggtree(nwk) %<+% nwkmeta +
  geom_tippoint(
    aes(
      color = country_of_origin,
      shape = serovar
    ),
    size = 3,
    alpha = 0.7
  ) +
  geom_tiplab(size = 2)

g


countrymeta <- as.data.frame(nwkmeta[, c("assembly", "country_of_origin")])
countrymeta

rownames(countrymeta) <- countrymeta$assembly
#View(countrydf)

#attach metadata to the treefile
g <- ggtree(nwk) %<+% nwkmeta + #here we are attaching the reordered metadata to the treefile
  geom_tippoint(aes(color=country_of_origin, shape=serovar), size=3, alpha=0.7) +
  geom_tiplab(size=2)

g

#these daylight trees might take a few minutes for R to compile
ggtree(nwk) %<+% nwkmeta + #here we are attaching the reordered metadata to the treefile
  geom_tippoint(aes(color=as.factor(mlst), shape=source_type_2), size=3, alpha=0.7) + 
  geom_tiplab(aes(label=assembly) , size=4) 

nwkmeta
source_shapes <- c(
  "Dairy" = 15,
  "Environmental" = 3,
  "Human Clinical" = 8,
  "Ovis aries" = 17,
  "Oyster" = 19,
  "Raw Animal Feed" = 12,
  "Unknown Food" = 25
)

ggtree(nwk) %<+% nwkmeta + #here we are attaching the reordered metadata to the treefile
  geom_tippoint(aes(color=serovar, shape=source_type_2), size=4, 
                alpha=0.8, #to change transparency of the shapes
                stroke = 2 #to increase the line-width of the points
                ) + 
  geom_tiplab(aes(label=assembly) , size=4) +
  scale_shape_manual(values = source_shapes) +
  scale_color_brewer(palette = "Set1")


#do we need to rename the tree tip labels maybe?
head(nwkmeta)
renamedtree <- nwk
renamedtree

#if we rename the tree, we need to reorder the metadata so the first column is the new labels
#column 9 is our new tree-tip label, so put it first. the order doesnt matter for the rest
nwkmeta2 <- nwkmeta

serovarcolors <- c("red","blue","green")

#plot the tree with the tip shapes as the source, and the tip color as the serovar
p <- ggtree(renamedtree) %<+% nwkmeta2 + #here we are attaching the reordered metadata to the treefile
  geom_tippoint(aes(color=serovar, shape=source_type_2), size=3, 
                alpha=0.8, #to change transparency of the shapes
                stroke = 2 #to increase the line-width of the points
  ) + 
  geom_tiplab(size=4,
              align = TRUE #use this to align the tip labels
              ) +
  scale_shape_manual(values = source_shapes) +
  scale_color_brewer(palette = "Set1") +
  #scale_color_manual(values = serovarcolors) +
  theme_tree2() + #tree scale bar
  theme(legend.position = c(0.3,0.8), legend.background = element_blank()) +
  guides(
    fill = guide_legend(
      title = "Country of Origin",
      order = 1
    ),
    color = guide_legend(
      title = "Serotype",
      ncol = 3,
      order = 2
    ),
    shape = guide_legend(
      title = "Isolation Source",
      ncol = 2,
      order = 3
    )
  )

p

#now lets add a heatmap
#only need a subset of the data
head(nwkmeta2)

countryhmp <- as.data.frame(nwkmeta2[, c('assembly', 'country_of_origin')])
countryhmp
row.names(countryhmp) <- countryhmp$assembly
countryhmp <- countryhmp[, 'country_of_origin', drop = FALSE]
countryhmp

colnames(countryhmp) <- "Country"

p2 <- p + new_scale_fill()

country_colors <- c(
  "Belgium" = "sienna1",    # Using Red as primary
  "Canada" = "red2",
  "Denmark" = "lightskyblue",
  "France" = "purple3",
  "Ireland" = "darkolivegreen2",
  "Netherlands" = "#FDDA24",
  "New Zealand" = "green4",
  "Norway" = "orchid",
  "Switzerland" = "burlywood4",
  "United Kingdom" = "slategray3"
)


p2

#all of the legends moved over to the right but whatever
p3 <- gheatmap(p2, countryhmp, 
         width = 0.1, 
         offset = 1600, #this moves the heatmap farther/closer to the tree
         color = "black",#this is to get the box outlines for the heatmap
         colnames = TRUE) +
  scale_fill_manual(values = country_colors, 
                    name = "Country of Origin", 
                    guide = guide_legend(order = 1),
                    labels = function(x) 
                      ifelse(x == "Canada", "Canada*", x)
                    )


p3

#to save
dpi = 300 #this will affect the resolution of your image
jpeg(filename="/Users/noorshubair/Desktop/R_code_files/chapter4/Linear_tree_w_country.jpeg", width=10*dpi, height=12.5*dpi, res=dpi)
p3
dev.off()

#if you want to add any more heatmaps:


# Adding mlst heatmap -----------------------------------------------------

mlsthmp <- as.data.frame(nwkmeta2[, c('assembly', 'mlst')])
mlsthmp
row.names(mlsthmp) <- mlsthmp$assembly
mlsthmp <- mlsthmp[, 'mlst', drop = FALSE]
mlsthmp$mlst <- as.factor(mlsthmp$mlst)

p4 <- p3 + new_scale_fill()

p5 <- gheatmap(p4, mlsthmp, 
               width = 0.1, 
               offset = 2000, #this moves the heatmap farther/closer to the tree
               color = "black",#this is to get the box outlines for the heatmap
               colnames = TRUE) +
  scale_fill_brewer(palette = "Set3",name = "MLST")
  #scale_fill_manual(values = country_colors, name = "Country")

p5


# unrooted tree -----------------------------------------------------------

unrooted <- ggtree(nwk, layout="daylight") %<+% nwkmeta + #here we are attaching the reordered metadata to the treefile
  geom_tippoint(aes(color=serovar, shape=source_type_2,
                    #put the x= and y= command here to spread the tip points
                    ), size=4, 
                alpha=0.8, #to change transparency of the shapes
                stroke = 2 #to increase the line-width of the points
  ) + 
  #geom_tiplab(aes(label=assembly) , size=4) +
  scale_shape_manual(name = "Isolation Source",
                     values = source_shapes) +
  scale_color_brewer(palette = "Set1") 


#to save
dpi = 300 #this will affect the resolution of your image
jpeg(filename="/Users/noorshubair/Desktop/R_code_files/chapter4/unrooted_tree.jpeg", width=10*dpi, height=12.5*dpi, res=dpi)
unrooted
dev.off()



#circular?
cp <- ggtree(renamedtree, layout = "circular", branch.length = "none") %<+% nwkmeta2 + #here we are attaching the reordered metadata to the treefile
  geom_tippoint(aes(color=serovar, shape=source_type_2), size=4, 
                alpha=0.8, #to change transparency of the shapes
                stroke = 2 #to increase the line-width of the points
  ) + 
  geom_tiplab2(size=4,
              #align = TRUE #use this to align the tip labels
  ) +
  scale_shape_manual(values = source_shapes) +
  scale_color_brewer(palette = "Set1")

cp2 <- cp + new_scale_fill()

gheatmap(cp2, countryhmp, 
         width = 0.1, offset = 7,
         color = "black",#this is to get the box outlines for the heatmap
         colnames = TRUE) +
  scale_fill_manual(values = country_colors)

