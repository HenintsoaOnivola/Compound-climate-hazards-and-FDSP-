library(viridis)
library(pacman)
library (DescTools)
library(raster)
library(rgdal)
library(tidyverse)
library(tibble)
library(ggplot2)
library(sf)
library(scales)
library(RColorBrewer)
library(cartography)
library (tidyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

HoA <- st_read("./qgis_files/input/world_shape2.geojson")

theme_MAP <-   theme_bw() +
  theme(title = element_text(size = 15),
        plot.title=element_text(hjust=0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size=10,face='bold'), 
        
  )

overlayRast <- function(rst) {
  r <- raster::crop(rst, HoA)
  rst <- raster::mask(r, HoA)
  tbl <- raster::rasterToPoints(rst, spatial = FALSE)
  tbl <- as_tibble(tbl)
  names(tbl) <- c('x', 'y', 'value')
  tbl
}


formatRast <- function(rst) {
  r_rotated <- raster::rotate(rst)
  r_croped <- raster::crop(r_rotated, HoA)
  rst <- raster::mask(r_croped, HoA)
  rst
}

Visu_raster<-function (output_name,the_min,the_max,title,bar_title,color_palette,all_ticks,all_labels,the_tbl){
  all_val<-na.omit(sort(unique(the_tbl$value)))
  my_cols<-color_palette
  if (all(all_val >= -0.01)){
    res<-seq(0,1,1/length(my_cols))
  }else
  {
    my_cols<- c('#767676','#FFFFFF',my_cols)
    positive_val=all_val[all_val>0]
    val_after_zero=positive_val[1]
    scaling_colormap<-c(the_min,0,seq(val_after_zero,the_max,(the_max-val_after_zero)/(length(my_cols)-3)))
    res<-rescale(scaling_colormap)
  }
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  print(ggplot() +
          geom_tile(data = the_tbl, aes(x = x, y = y, fill = value)) +
          geom_sf(data = HoA, fill = NA,lwd=0.5) +
          scale_fill_gradientn(colors = my_cols,name=bar_title, limits=c(the_min,the_max),values=res,breaks=all_ticks,labels=all_labels)+   
          labs(title = title) +
          theme_MAP)
  dev.off()
}

Visu_raster_discrete <-function (output_name,title,bar_title,color_palette,the_lim,all_ticks,all_labels,the_tbl){
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  print(ggplot() +
          geom_tile(data = the_tbl, aes(x = x, y = y, fill = factor(value))) +
          geom_sf(data = HoA, fill = NA,lwd=0.5) +
          scale_fill_manual(values=color_palette, name=bar_title,limits=the_lim,breaks=all_ticks,labels=all_labels) +   
          labs(title = title) +
          theme_MAP)
  dev.off()
}

dominant_hazard <- function(the_stack) {
  ind_max<-which.max(the_stack)
  if (identical(ind_max,integer(0))){
    ind_max=NA
  }
  ind_max
}

##transforming raw indices to be comparable (normalization and distribution shift)
#######################################################################
raw_heat_hist<-raster("./raw_indices/da_n_day_hist.tif")
raw_heat_ssp<-raster("./raw_indices/da_n_day_ssp.tif")
raw_drought_hist<-raster("./raw_indices/SPEI_low10_hist.tif")
raw_drought_ssp<-raster("./raw_indices/SPEI_low10_ssp.tif")
raw_flood_hist<-raster("./raw_indices/SDII_hist.tif")
raw_flood_ssp<-raster("./raw_indices/SDII_ssp.tif")

raw_heat_hist<-formatRast(raw_heat_hist)
raw_heat_ssp<-formatRast(raw_heat_ssp)
raw_drought_hist<-formatRast(raw_drought_hist)
raw_drought_ssp<-formatRast(raw_drought_ssp)
raw_flood_hist<-formatRast(raw_flood_hist)
raw_flood_ssp<-formatRast(raw_flood_ssp)

##baseline
raw_heat_hist2<-raw_heat_hist
raw_heat_hist2[raw_heat_hist2<1]<-NA
heat_hist_nrm<-(raw_heat_hist2-minValue(raw_heat_hist2))/(maxValue(raw_heat_hist2)-minValue(raw_heat_hist2))
raw_drought_hist2<- -raw_drought_hist
drought_hist_nrm<- (raw_drought_hist2-minValue(raw_drought_hist2))/(maxValue(raw_drought_hist2)-minValue(raw_drought_hist2))
flood_hist_nrm<-(raw_flood_hist-minValue(raw_flood_hist))/(maxValue(raw_flood_hist)-minValue(raw_flood_hist))
##future
raw_heat_ssp2<-raw_heat_ssp
raw_heat_ssp2[raw_heat_ssp2<1]<-NA
heat_ssp_nrm<-(raw_heat_ssp2-minValue(raw_heat_ssp2))/(maxValue(raw_heat_ssp2)-minValue(raw_heat_ssp2))
raw_drought_ssp2<- -raw_drought_ssp
drought_ssp_nrm<- (raw_drought_ssp2-minValue(raw_drought_ssp2))/(maxValue(raw_drought_ssp2)-minValue(raw_drought_ssp2))
flood_ssp_nrm<-(raw_flood_ssp-minValue(raw_flood_ssp))/(maxValue(raw_flood_ssp)-minValue(raw_flood_ssp))

#plotting the distribution of raw values
tiff('./summary_figures_bis/raw_hazard_distribution_hist.tiff', units="in", width=7, height=4, res=500,compression = 'lzw')
par(mfrow = c(2, 3))
hist(raw_heat_hist2,xlab='number if days with HI>41°C',ylab='frequency',breaks = "Scott",main='')
title(main = "heat",adj = 0)
abline(v=quantile(raw_heat_hist2, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=30,col='red',lwd=1.5,lty=2)
hist(raw_drought_hist2,xlab='-(SPEI)',ylab='frequency',breaks=75,main='')
title(main = "drought",adj = 0)
abline(v=quantile(raw_drought_hist2, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=1,col='red',lwd=1.5,lty=2)
mtext("Baseline", side = 3, line = 3, outer =FALSE)
hist(raw_flood_hist,xlab='SDII',ylab='frequency',breaks = "Scott",main='')
title(main = "flood",adj = 0)
abline(v=quantile(raw_flood_hist, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=10,col='red',lwd=1.5,lty=2)
legend('topright',inset=c(-0.1, -0.5), legend=c("80th percentile", "severe threshold"),
       col=c("blue", "red"), lty=1:2, cex=0.8,xpd=TRUE)

hist(raw_heat_ssp2,xlab='number if days with HI>41°C',ylab='frequency',breaks = "Scott",main='')
title(main = "heat",adj = 0)
abline(v=quantile(raw_heat_ssp2, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=30,col='red',lwd=1.5,lty=2)
hist(raw_drought_ssp2,xlab='-(SPEI)',ylab='frequency',breaks=80,main='')
title(main = "drought",adj = 0)
abline(v=quantile(raw_drought_ssp2, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=1,col='red',lwd=1.5,lty=2)
mtext("Future", side = 3, line = 3, outer = FALSE)
hist(raw_flood_ssp,xlab='SDII',ylab='frequency',breaks = "Scott",main='')
title(main = "flood",adj = 0)
abline(v=quantile(raw_flood_ssp, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=10,col='red',lwd=1.5,lty=2)

dev.off()

#set the threshold to be considered as 'extreme' and 'severe' (class 4  and class 5 hazard), and normalize their values
heat_thresh_nrm_hist<-(30-minValue(raw_heat_hist2))/(maxValue(raw_heat_hist2)-minValue(raw_heat_hist2))
drought_thresh_nrm_hist<- (1-minValue(raw_drought_hist2))/(maxValue(raw_drought_hist2)-minValue(raw_drought_hist2))
flood_thresh_nrm_hist<- (10-minValue(raw_flood_hist))/(maxValue(raw_flood_hist)-minValue(raw_flood_hist))
heat_thresh_nrm_ssp<-(30-minValue(raw_heat_ssp))/(maxValue(raw_heat_ssp)-minValue(raw_heat_ssp))
drought_thresh_nrm_ssp<- (1-minValue(raw_drought_ssp2))/(maxValue(raw_drought_ssp2)-minValue(raw_drought_ssp2))
flood_thresh_nrm_ssp<- (10-minValue(raw_flood_ssp))/(maxValue(raw_flood_ssp)-minValue(raw_flood_ssp))
#plotting distribution frequency of the hazards
theme_frequency<- theme(title = element_text(size = 15),
                        plot.title=element_text(hjust=0.5),
                        axis.text=element_text(size=10),
                        axis.line = element_line(colour = "black"),
                        plot.background = element_blank(),
                        panel.background=element_rect(fill = "white"),
                        panel.border = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.title = element_text(size=10,face='bold'),
                        legend.key = element_rect(fill = "white"))
df_freq_hist<- data.frame(x = values(heat_hist_nrm), y = values(drought_hist_nrm),z=values(flood_hist_nrm)) %>%
  gather(key, value)
tiff('./summary_figures_bis/hazard_distribution_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_hist, aes(value, y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = heat_thresh_nrm_hist, linetype="dashed", 
             color = "red")+
  geom_vline(xintercept = drought_thresh_nrm_hist, linetype="dashed", 
             color = "yellow")+
  geom_vline(xintercept = flood_thresh_nrm_hist, linetype="dashed", 
             color = "blue")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - baseline',x='hazard index value',y='density')
dev.off()

df_freq_ssp<- data.frame(x = values(heat_ssp_nrm), y = values(drought_ssp_nrm),z=values(flood_ssp_nrm)) %>%
  gather(key, value)
tiff('./summary_figures_bis/hazard_distribution_ssp.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_ssp, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = heat_thresh_nrm_ssp, linetype="dashed", 
             color = "red")+
  geom_vline(xintercept = drought_thresh_nrm_ssp, linetype="dashed", 
             color = "yellow")+
  geom_vline(xintercept = flood_thresh_nrm_ssp, linetype="dashed", 
             color = "blue")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - future',x='hazard index value',y='density')
dev.off()

##shifting the hazard distributions so that the threshold for 'severe' coincide to 0.6 (high and extremely high)
#baseline
heat_hist_sft<-heat_hist_nrm+(0.6-heat_thresh_nrm_hist)
drought_hist_sft<-drought_hist_nrm+(0.6-drought_thresh_nrm_hist)
flood_hist_sft<-flood_hist_nrm+(0.6-flood_thresh_nrm_hist)
heat_hist_sft[heat_hist_sft>=1]<-1
drought_hist_sft[drought_hist_sft>=1]<-1
flood_hist_sft[flood_hist_sft>=1]<-1
heat_hist_sft[heat_hist_sft<=0]<-0
drought_hist_sft[drought_hist_sft<=0]<-0
flood_hist_sft[flood_hist_sft<=0]<-0

heat_thresh_sft_hist<- heat_thresh_nrm_hist + (0.6-heat_thresh_nrm_hist) 
drought_thresh_sft_hist<- drought_thresh_nrm_hist + (0.6-drought_thresh_nrm_hist) 
flood_thresh_sft_hist<- flood_thresh_nrm_hist + (0.6-flood_thresh_nrm_hist)

#future
heat_ssp_sft<-heat_ssp_nrm+(0.6-heat_thresh_nrm_ssp)
drought_ssp_sft<-drought_ssp_nrm+(0.6-drought_thresh_nrm_ssp)
flood_ssp_sft<-flood_ssp_nrm+(0.6-flood_thresh_nrm_ssp)
heat_ssp_sft[heat_ssp_sft>=1]<-1
drought_ssp_sft[drought_ssp_sft>=1]<-1
flood_ssp_sft[flood_ssp_sft>=1]<-1
heat_ssp_sft[heat_ssp_sft<=0]<-0
drought_ssp_sft[drought_ssp_sft<=0]<-0
flood_ssp_sft[flood_ssp_sft<=0]<-0

heat_thresh_sft_ssp<- heat_thresh_nrm_ssp + (0.6-heat_thresh_nrm_ssp) 
drought_thresh_sft_ssp<- drought_thresh_nrm_ssp + (0.6-drought_thresh_nrm_ssp) 
flood_thresh_sft_ssp<- flood_thresh_nrm_ssp + (0.6-flood_thresh_nrm_ssp)

#plot the new shifted distributions
df_freq_hist_sft<- data.frame(x = values(heat_hist_sft), y = values(drought_hist_sft),z=values(flood_hist_sft)) %>%
  gather(key, value)
tiff('./summary_figures_bis/hazard_distribution_shifted_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_hist_sft, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = heat_thresh_sft_hist, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = drought_thresh_sft_hist, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = flood_thresh_sft_hist, linetype="dashed", 
             color = "black")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - baseline',x='hazard index value',y='density')
dev.off()

df_freq_ssp_sft<- data.frame(x = values(heat_ssp_sft), y = values(drought_ssp_sft),z=values(flood_ssp_sft)) %>%
  gather(key, value)
tiff('./summary_figures_bis/hazard_distribution_shifted_ssp.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_ssp_sft, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = heat_thresh_sft_ssp, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = drought_thresh_sft_ssp, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = flood_thresh_sft_ssp, linetype="dashed", 
             color = "black")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - future',x='hazard index value',y='density')
dev.off()

heat_hist_index<-heat_hist_sft
heat_hist_index[raw_heat_hist<1]<-minValue(heat_hist_index)
drought_hist_index<-drought_hist_sft
flood_hist_index<-flood_hist_sft

heat_ssp_index<-heat_ssp_sft
heat_ssp_index[raw_heat_ssp<1]<-minValue(heat_ssp_index)
drought_ssp_index<-drought_ssp_sft
flood_ssp_index<-flood_ssp_sft

heat_change_index<-heat_ssp_index-heat_hist_index
drought_change_index<-drought_ssp_index-drought_hist_index
flood_change_index<-flood_ssp_index-flood_hist_index


#write the indices into raster files
writeRaster(heat_hist_index, filename = ("./qgis_files/output_bis/heat_hist_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(drought_hist_index, filename = ("./qgis_files/output_bis/drought_hist_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(flood_hist_index, filename = ("./qgis_files/output_bis/flood_hist_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(heat_ssp_index, filename = ("./qgis_files/output_bis/heat_ssp_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(drought_ssp_index, filename = ("./qgis_files/output_bis/drought_ssp_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(flood_ssp_index, filename = ("./qgis_files/output_bis/flood_ssp_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(heat_change_index, filename = ("./qgis_files/output_bis/heat_change_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(drought_change_index, filename = ("./qgis_files/output_bis/drought_change_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(flood_change_index, filename = ("./qgis_files/output_bis/flood_change_nrm.tif"), format = "GTiff",overwrite=TRUE)


##visualizing individual climate hazard indices (current, future, change)
##visualizing heat hazard
heat_hist_tbl <- overlayRast(heat_hist_index)
heat_ssp_tbl <- overlayRast(heat_ssp_index)
heat_change_index<-heat_ssp_index-heat_hist_index
heat_change_tbl<-overlayRast(heat_change_index)
the_min<-min(c(minValue(heat_hist_index),minValue(heat_ssp_index)))
the_max<-max(c(maxValue(heat_hist_index),maxValue(heat_ssp_index)))
col_map<-rev(viridis::rocket(27))
col_map<-col_map[4:(length(col_map)-1)]
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster("./summary_figures_bis/heat_hist_index.tiff",color_bar_min,color_bar_max,"heat hazard - baseline",'heat index',col_map,the_ticks,the_labels,heat_hist_tbl)
Visu_raster("./summary_figures_bis/heat_ssp_index.tiff",color_bar_min,color_bar_max,"heat hazard - future",'heat index',col_map,the_ticks,the_labels,heat_ssp_tbl)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster("./summary_figures_bis/heat_change_index.tiff",minValue(heat_change_index),maxValue(heat_change_index),"heat hazard - change",'heat change',col_map,waiver(),waiver(),heat_change_tbl)

##visualizing drought hazard
drought_hist_tbl <- overlayRast(drought_hist_index)
drought_ssp_tbl <- overlayRast(drought_ssp_index)
drought_change_index<-drought_ssp_index-drought_hist_index
drought_change_tbl<-overlayRast(drought_change_index)
the_min<-min(c(minValue(drought_hist_index),minValue(drought_ssp_index)))
the_max<-max(c(maxValue(drought_hist_index),maxValue(drought_ssp_index)))
col_map<-carto.pal('sand.pal',20)
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster("./summary_figures_bis/drought_hist_index.tiff",color_bar_min,color_bar_max,"drought hazard - baseline",'drought index',col_map,the_ticks,the_labels,drought_hist_tbl)
Visu_raster("./summary_figures_bis/drought_ssp_index.tiff",color_bar_min,color_bar_max,"drought hazard - future",'drought index',col_map,the_ticks,the_labels,drought_ssp_tbl)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster("./summary_figures_bis/drought_change_index.tiff",minValue(drought_change_index),maxValue(drought_change_index),"drought hazard - change",'drought change',col_map,waiver(),waiver(),drought_change_tbl)

##visualizing flood hazard
flood_hist_tbl <- overlayRast(flood_hist_index)
flood_ssp_tbl <- overlayRast(flood_ssp_index)
flood_change_index<-flood_ssp_index-flood_hist_index
flood_change_tbl<-overlayRast(flood_change_index)
the_min<-min(c(minValue(flood_hist_index),minValue(flood_ssp_index)))
the_max<-max(c(maxValue(flood_hist_index),maxValue(flood_ssp_index)))
col_map<-carto.pal('blue.pal',20)
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster("./summary_figures_bis/flood_hist_index.tiff",color_bar_min,color_bar_max,"flood hazard - baseline",'flood index',col_map,the_ticks,the_labels,flood_hist_tbl)
Visu_raster("./summary_figures_bis/flood_ssp_index.tiff",color_bar_min,color_bar_max,"flood hazard - future",'flood index',col_map,the_ticks,the_labels,flood_ssp_tbl)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster("./summary_figures_bis/flood_change_index.tiff",minValue(flood_change_index),maxValue(flood_change_index),"flood hazard - change",'flood change',col_map,waiver(),waiver(),flood_change_tbl)

##############################################################
##plotting locations of 'severe' hazard

high_heat_hist<-(heat_hist_index>=0.6)
high_drought_hist<-(drought_hist_index>=0.6)
high_flood_hist<-(flood_hist_index>=0.6)

high_heat_hist_tbl<-overlayRast(high_heat_hist)
high_drought_hist_tbl<-overlayRast(high_drought_hist)
high_flood_hist_tbl<-overlayRast(high_flood_hist)

high_heat_ssp<-(heat_ssp_index>=0.6)
high_drought_ssp<-(drought_ssp_index>=0.6)
high_flood_ssp<-(flood_ssp_index>=0.6)

high_heat_ssp_tbl<-overlayRast(high_heat_ssp)
high_drought_ssp_tbl<-overlayRast(high_drought_ssp)
high_flood_ssp_tbl<-overlayRast(high_flood_ssp)

my_pals=c('white','black')
Visu_raster_discrete("./summary_figures_bis/high_heat_hist.tiff","severe heat - baseline",'heat stress',my_pals,c('0','1'),c(0,1),c('','high'),high_heat_hist_tbl)
Visu_raster_discrete("./summary_figures_bis/high_drought_hist.tiff","severe drought- baseline",'drought intensity',my_pals,c('0','1'),c(0,1),c('','high'),high_drought_hist_tbl)
Visu_raster_discrete("./summary_figures_bis/high_flood_hist.tiff","severe flooding - baseline",'flood intensity',my_pals,c('0','1'),c(0,1),c('','high'),high_flood_hist_tbl)
Visu_raster_discrete("./summary_figures_bis/high_heat_ssp.tiff","severe heat - future",'heat stress',my_pals,c('0','1'),c(0,1),c('','high'),high_heat_ssp_tbl)
Visu_raster_discrete("./summary_figures_bis/high_drought_ssp.tiff","severe drought- future",'drought intensity',my_pals,c('0','1'),c(0,1),c('','high'),high_drought_ssp_tbl)
Visu_raster_discrete("./summary_figures_bis/high_flood_ssp.tiff","severe flooding - future",'flood intensity',my_pals,c('0','1'),c(0,1),c('','high'),high_flood_ssp_tbl)

##showing all individual severe hazards with their combinations
high_ind_hist <- high_heat_hist
high_ind_hist[,]<-0
high_ind_hist[high_heat_hist==1 & high_drought_hist==0 & high_flood_hist==0]<-1
high_ind_hist[high_heat_hist==0 & high_drought_hist==1 & high_flood_hist==0]<-2
high_ind_hist[high_heat_hist==0 & high_drought_hist==0 & high_flood_hist==1]<-3
high_ind_hist[high_heat_hist==1 & high_drought_hist==1 & high_flood_hist==0]<-4
high_ind_hist[high_heat_hist==1 & high_drought_hist==0 & high_flood_hist==1]<-5
high_ind_hist[high_heat_hist==0 & high_drought_hist==1 & high_flood_hist==1]<-6
high_ind_hist[high_heat_hist==1 & high_drought_hist==1 & high_flood_hist==1]<-7
high_ind_hist_tbl<-overlayRast(high_ind_hist)

high_ind_ssp <- high_heat_ssp
high_ind_ssp[,]<-0
high_ind_ssp[high_heat_ssp==1 & high_drought_ssp==0 & high_flood_ssp==0]<-1
high_ind_ssp[high_heat_ssp==0 & high_drought_ssp==1 & high_flood_ssp==0]<-2
high_ind_ssp[high_heat_ssp==0 & high_drought_ssp==0 & high_flood_ssp==1]<-3
high_ind_ssp[high_heat_ssp==1 & high_drought_ssp==1 & high_flood_ssp==0]<-4
high_ind_ssp[high_heat_ssp==1 & high_drought_ssp==0 & high_flood_ssp==1]<-5
high_ind_ssp[high_heat_ssp==0 & high_drought_ssp==1 & high_flood_ssp==1]<-6
high_ind_ssp[high_heat_ssp==1 & high_drought_ssp==1 & high_flood_ssp==1]<-7
high_ind_ssp_tbl<-overlayRast(high_ind_ssp)

my_pals=c('0'='white','1'='orangered1','2'='goldenrod3','3'='steelblue3','4'='saddlebrown','5'='olivedrab1','6'='purple','7'='grey')

Visu_raster_discrete("./summary_figures_bis/high_ind_hist.tiff","severe hazard - baseline",'hazard type',my_pals,c('0','1','2','3','4','5','6','7'),c('0','1','2','3','4','5','6','7'),c('','heat','drought','flood','heat and drougth','heat and flood','drought and flood','all'),high_ind_hist_tbl)
Visu_raster_discrete("./summary_figures_bis/high_ind_ssp.tiff","severe hazard - future",'hazard type',my_pals,c('0','1','2','3','4','5','6','7'),c('0','1','2','3','4','5','6','7'),c('','heat','drought','flood','heat and drougth','heat and flood','drought and flood','all'),high_ind_ssp_tbl)

###Plotting hazard per region
regional_files<-c('AMERICAs_15m_UNHCR.gpkg',"ASIA_15m_UNHCR.gpkg","East_and_horn_of_Africa_shapefile.gpkg",
                 "Europe_new.gpkg","MENA_polygon_2.gpkg","Southern_Africa_Ash.gpkg","West_and_central_new.gpkg")


all_region<-c('Americas','Asia','East and Horn of Africa','Europe','MENA','Southern Africa','West and Central Africa')
all_types=c('heat','drought','flood','heat and drougth','heat and flood','drought and flood','all')
all_colors=c('orangered1','goldenrod3','steelblue3','saddlebrown','olivedrab1','purple','grey')
##baseline
percentage=c()
type=c()
region=c()
for (count_region in 1:length(regional_files)){
  polygon_rst <- st_read(paste0("./qgis_files/input/",regional_files[count_region]))
  r <- raster::crop(high_ind_hist, polygon_rst)
  rst_ini <- raster::mask(r, polygon_rst)
  rst_all_values<-rst_ini
  rst_all_values[rst_ini>=0]=1
  rst_all_high<-rst_ini
  rst_all_high[rst_ini>=1]=1
  all_type_high<-cellStats(rst_all_high, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
  for (type_ind in (1:length(all_types))){
    type_rast<-rst_ini
    type_rast[rst_ini!=type_ind]=0
    type_rast[type_rast==type_ind]=1
    proportion_type<-cellStats(type_rast, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
    
    percentage<-c(percentage,proportion_type)
    type<-c(type,all_types[type_ind])
    region<-c(region,all_region[count_region])
  }
}
all_colors=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3','heat and drougth'='saddlebrown','heat and flood'='olivedrab1','drought and flood'='purple','all'='grey')
the_df<-data.frame(percentage,type,region)
the_df$region<-factor(the_df$region)
the_df2<-the_df[the_df$percentage!=0,]
tiff('./summary_figures_bis/hazard_per_region_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(the_df,aes(x=region,y=percentage,fill=type))+
  geom_bar(stat='identity')+
  labs(y='Percentage of area under severe hazard')+
  #scale_fill_manual(values=all_colors,breaks=all_types,labels=all_types)+
  scale_fill_manual(values=all_colors,breaks=unique(the_df2$type),labels=all_types[match(unique(the_df2$type),all_types)])+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  ylim(0,1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  ggtitle('severe hazard - baseline')
dev.off()

##future
percentage=c()
type=c()
region=c()
for (count_region in 1:length(regional_files)){
  polygon_rst <- st_read(paste0("./qgis_files/input/",regional_files[count_region]))
  r <- raster::crop(high_ind_ssp, polygon_rst)
  rst_ini <- raster::mask(r, polygon_rst)
  rst_all_values<-rst_ini
  rst_all_values[rst_ini>=0]=1
  rst_all_high<-rst_ini
  rst_all_high[rst_ini>=1]=1
  all_type_high<-cellStats(rst_all_high, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
  for (type_ind in (1:length(all_types))){
    type_rast<-rst_ini
    type_rast[rst_ini!=type_ind]=0
    type_rast[type_rast==type_ind]=1
    proportion_type<-cellStats(type_rast, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
    
    percentage<-c(percentage,proportion_type)
    type<-c(type,all_types[type_ind])
    region<-c(region,all_region[count_region])
  }
}
all_colors=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3','heat and drougth'='saddlebrown','heat and flood'='olivedrab1','drought and flood'='purple','all'='grey')
the_df<-data.frame(percentage,type,region)
the_df$region<-factor(the_df$region)
the_df2<-the_df[the_df$percentage!=0,]
tiff('./summary_figures_bis/hazard_per_region_ssp.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(the_df,aes(x=region,y=percentage,fill=type))+
  geom_bar(stat='identity')+
  labs(y='Percentage of area under severe hazard')+
  #scale_fill_manual(values=all_colors,breaks=all_types,labels=all_types)+
  scale_fill_manual(values=all_colors,breaks=unique(the_df2$type),labels=all_types[match(unique(the_df2$type),all_types)])+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  ylim(0,1)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  ggtitle('severe hazard - future')
dev.off()


#dominant hazard per point locations
my_pals=c('orangered1','goldenrod3','steelblue3')
my_labels=c('heat','drought','flood')
#baseline
hist_stack<-stack(heat_hist_index,drought_hist_index,flood_hist_index)
max_ind_hist <- calc(hist_stack,dominant_hazard)
max_ind_hist_tbl<-overlayRast(max_ind_hist)
Visu_raster_discrete('./summary_figures_bis/common_hazard_hist.tiff','most common hazard - baseline','hazard type',my_pals,c('1','2','3'),c('1','2','3'),my_labels,max_ind_hist_tbl)
#future
ssp_stack<-stack(heat_ssp_index,drought_ssp_index,flood_ssp_index)
max_ind_ssp<- calc(ssp_stack,dominant_hazard)
max_ind_ssp_tbl<-overlayRast(max_ind_ssp)
Visu_raster_discrete('./summary_figures_bis/common_hazard_ssp.tiff','most common hazard - future','hazard type',my_pals,c('1','2','3'),c('1','2','3'),my_labels,max_ind_ssp_tbl)




