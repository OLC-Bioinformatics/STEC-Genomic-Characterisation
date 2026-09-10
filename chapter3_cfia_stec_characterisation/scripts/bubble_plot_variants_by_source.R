#bubble plots to show percentage from source for each profile/variant

#load required packages
library(ggtree)
library(reshape2)
library(ggplot2)
library(ggimage)
library(readxl)
library(RColorBrewer)
library(ggsci)

#set our working directory
setwd("/Users/noorshubair/Desktop/R_code_files/chapter3")

#read in our excel data
chpt3stxvariantsdata <- read_excel("CFIA_STEC.xlsx", sheet="R_bubbleplot")
chpt3stxvariantsdata

#pull the source columns into one column with Source and total_from_source_with_variant as headers
library(tidyr)
chpt3stxvariantsdata2 <- pivot_longer(
  data = chpt3stxvariantsdata,
  cols = c(3:9),
  names_to = "Source",
  values_to = "total_from_source_with_variant"
)
chpt3stxvariantsdata2

#theme with nice horizontal dots along axis
theme_dotplotx <- function() {
  theme( ## remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    panel.grid.minor.x = element_blank() ,
    ## explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(color="grey70", linetype = 3),
    axis.text.y = element_text(size=rel(1.2)),
    ## use a white backgrounsd
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"))
}

#a custom colour palette
source_colours <- c("yellow","orange","white",#Dairy,Env, Flour
                    #"pink",#lamb
                    "red","purple1",#meat,nuts/seed
                    "grey50","green3" #other, vegetable
)

#add percentage info for 
library(dplyr)
chpt3stxvariantsdata2
chpt3stxvariantsdata2$percent <-  chpt3stxvariantsdata2$total_from_source_with_variant/chpt3stxvariantsdata2$Total_w_StxVariant*100
chpt3stxvariantsdata2


ggplot(chpt3stxvariantsdata2,
       aes(x=Stx_variant,#x=reorder(Stx_variant, Total_w_StxVariant), #reordering based on the number of sequences, so largest to smallest number
           y=percent, 
           group=Source, 
           size=total_from_source_with_variant, 
           fill=Source, 
           color="black")#use colour to change the bubble outline color
) +
  geom_point(aes(y=percent,
                 x=Stx_variant#x=reorder(Stx_variant, Total_w_StxVariant)
                 ), 
             stat="identity", alpha=0.5, shape=21, color="black") +
  scale_size(range=c(3,20), #use to change the range of bubble sizes
             name="No. Sequences",
             breaks=c(1,10,100)
  ) +
  scale_fill_manual(values = source_colours, #using the colour palette I created above
                    #to manually change the labels in our legend:
                    #labels = c("Animal","Animal Feed","Environmental","Food",
                    #           "Human Clinical","Not Provided","Other")
                    ) +
  scale_y_continuous(limits = c(0,100)) +
  theme_dotplotx() +
  #now use below to change the text sizes for the different variables
  theme(axis.text.x=element_text(size = 15), axis.title.x = element_text(size = 17),
        axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 17),
        legend.text = element_text(size = 12),legend.title = element_text(size = 15)) +
  labs(x="Stx Variant", y="Proportion from Source (%)") +
  coord_flip() + #swapping the x and y axes
  #legend fixes
  #use guides below to change the size of the points in the legend
  guides(fill = guide_legend(override.aes = list(size=5)), 
         size=guide_legend(override.aes = list(stroke=1)))


#ordered from most common to least
bubblesource <- ggplot(chpt3stxvariantsdata2,
       aes(x=reorder(Stx_variant, Total_w_StxVariant), #reordering based on the number of sequences, so largest to smallest number
           y=percent, 
           group=Source, 
           size=total_from_source_with_variant, 
           fill=Source, 
           color="black")#use colour to change the bubble outline color
) +
  geom_point(aes(y=percent,
                 x=reorder(Stx_variant, Total_w_StxVariant)
  ), 
  stat="identity", alpha=0.7, shape=21, color="black") +
  scale_size(range=c(3,15), #use to change the range of bubble sizes
             name="No. Sequences",
             breaks=c(1,10,100)
  ) +
  scale_fill_manual(values = source_colours, #using the colour palette I created above
                    #to manually change the labels in our legend:
                    #labels = c("Animal","Animal Feed","Environmental","Food",
                    #           "Human Clinical","Not Provided","Other")
  ) +
  scale_y_continuous(limits = c(0,100)) +
  theme_dotplotx() +
  #now use below to change the text sizes for the different variables
  theme(axis.text.x=element_text(size = 15), axis.title.x = element_text(size = 17),
        axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 17),
        legend.text = element_text(size = 12),legend.title = element_text(size = 15)) +
  labs(x="Stx Variant", y="Proportion from Source (%)") +
  coord_flip() + #swapping the x and y axes
  #legend fixes
  #use guides below to change the size of the points in the legend
  guides(fill = guide_legend(override.aes = list(size=5)), 
         size=guide_legend(override.aes = list(stroke=1)))


dpi=300
jpeg("Variant_by_source.jpg", heigh=21*dpi, width=14*dpi, res=dpi)
bubblesource
dev.off()

