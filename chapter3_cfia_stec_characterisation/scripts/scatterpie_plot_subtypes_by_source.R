#scatterpie plot for stx protein carriage by serotype
#scatterpie plot for tellurite data
# Install devtools if not already installed
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}

# Install scatterpie from GitHub
#devtools::install_github("GuangchuangYu/scatterpie")

library(scatterpie)
library(readxl)

setwd("/Users/noorshubair/Desktop/R_code_files/chapter3")

#theme for plot
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

#stxsource_df <- read_excel("CFIA_STEC.xlsx", sheet="stx_source_pies_R")
stxsource_df <- read_excel("CFIA_STEC.xlsx", sheet="stx_pies_497")
stxsource_df


custom_x_labels <- c("Stx1a" = 5, 
                     "Stx2a"=35,
                     "Stx2d"=60,
                     "Stx2c"=85,
                     "Stx1c"=110,
                     "Stx2b"=135,
                     "Stx2k"=160,
                     "Stx2e"=185,
                     "Stx2g"=210,
                     "Stx2i"=235,
                     "Stx1d"=260,
                     "Stx2h"=285,
                     "Stx2j"=310
                     )

head(stxsource_df)

source_colours <- c("yellow","orange","white",#Dairy,Env, Flour
                    #"pink",#lamb
                    "red","purple1",#meat,nuts/seed
                    "grey50","green3" #other, vegetable
                    )

# larger pies -------------------------------------------------------------

stxsource_df <- read_excel("CFIA_STEC.xlsx", sheet="stx_pies_497")

custom_x_labels <- c("Stx1a" = 5, 
                     "Stx2a"=40,
                     "Stx2d"=75,
                     "Stx2c"=110,
                     "Stx1c"=145,
                     "Stx2b"=180,
                     "Stx2k"=215,
                     "Stx2e"=250,
                     "Stx2g"=285,
                     "Stx2i"=320,
                     "Stx1d"=355,
                     "Stx2h"=390,
                     "Stx2j"=425
)

source_colours <- c("yellow","orange","white",#Dairy,Env, Flour
                    #"pink",#lamb
                    "red","purple1",#meat,nuts/seed
                    "grey80","green3" #other, vegetable
)

plot <- ggplot() + 
  geom_scatterpie2(aes(x=x_axis2, y=y_axis, 
                       group=Source,
                       r=radius2 #the radius of the pies
  ), data=stxsource_df,
  alpha = 0.8,
  linewidth = 0.1,
  cols="Source",
  legend_name = "Source",
  long_format=TRUE,
  label_show_ratio = TRUE, 
  #label_radius = 1.1,
  #label_threshold = 1, fontsize = 1,
  donut_radius=0.3 #you HAVE to set a radius in your data table for this to work
  ) + 
  #facet_wrap(.~Facet, nrow=2) +
  coord_equal() +
  scale_x_continuous(labels = names(custom_x_labels), 
                     breaks = custom_x_labels) +
  scale_y_continuous(limits = c(-17,100),
                     breaks = c(1,25,50,75,100)) +
  scale_fill_manual(values = source_colours,
                    #labels = c("O157:H7","O26","O103","O111","O121","O145","Other")
  ) +
  #theme_bw() +
  theme_dotplotx() +
  theme(axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17),
        axis.text.x = element_text(size = 15, 
                                   #angle=60, 
                                   vjust=0.7),#, face = "italic"), 
        #axis.title.x = element_text(size = 17, vjust = 0.1),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 14), legend.title = element_text(size = 15),
        strip.background.x = element_rect(fill = "grey95", color="black"),
        strip.text.x = element_text(size = 15)
  ) +
  geom_text(data = stxsource_df, 
            aes(x=x_axis2, y=y_axis, label=total_protein),
            hjust=0.5, vjust=0.5, size = 4) +
  labs(x = "Protein", y="Genomes Encoding Stx Subtype (%)") +
  guides(fill = guide_legend(nrow = 1))

dpi=300
jpeg("stx_by_source_pies_PRJNA454819.jpeg",width=15*dpi, height=8*dpi, res=dpi)
plot
dev.off()


# old codes ---------------------------------------------------------------



#lets specify the order of the faceting variable so that the including O157 is first.
stxsource_df$Facet <- factor(stxsource_df$Facet, 
                             levels = c("Including O157", "Excluding O157"))
stxsource_df

pieplotstx <- 
  ggplot() + 
  geom_scatterpie2(aes(x=x_axis, y=y_axis, 
                       group=Source,
                       r=radius #the radius of the pies
  ), data=stxsource_df,
  alpha = 0.8,
  cols="Source",
  legend_name = "Source",
  long_format=TRUE,
  label_show_ratio = TRUE, 
  #label_radius = 1.1,
  #label_threshold = 1, fontsize = 1,
  donut_radius=0.3 #you HAVE to set a radius in your data table for this to work
  ) + 
  #facet_wrap(.~Facet, nrow=2) +
  coord_equal() +
  scale_x_continuous(labels = names(custom_x_labels), 
                     breaks = custom_x_labels) +
  scale_y_continuous(limits = c(-15,100),
                     breaks = c(1,25,50,75,100)) +
  scale_fill_manual(values = source_colours,
                    #labels = c("O157:H7","O26","O103","O111","O121","O145","Other")
                    ) +
  #theme_bw() +
  theme_dotplotx() +
  theme(axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13, angle=60, vjust=0.8),#, face = "italic"), 
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        strip.background.x = element_rect(fill = "grey95", color="black"),
        strip.text.x = element_text(size = 15)
        ) +
  geom_text(data = stxsource_df, 
            aes(x=x_axis, y=y_axis, label=total_protein),
            hjust=0.5, vjust=0.5, size = 4) +
    labs(x = "Protein", y="Genomes Encoding Protein (%)")

pieplotstx

#save the image
dpi=300
jpeg("stx_by_source_pies.jpeg",width=13*dpi, height=10*dpi, res=dpi)
pieplotstx
dev.off()




