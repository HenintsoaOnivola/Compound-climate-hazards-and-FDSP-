# This script transforms raw climate data into climate hazard indices
# and create visualization of the individual climate hazards distribution

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
library (colorspace)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#inpu the word shape file containing country boundaries
HoA <- st_read("./input_data/spatial_data/world_shape2.geojson")

##############################################
#defining all functions and plot formating
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
  # r_rotated <- raster::rotate(rst)  ## already rotated
  r_croped <- raster::crop(rst, HoA)
  rst2 <- raster::mask(r_croped, HoA)
  rst2
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
#####################################################

##transforming raw indices to be comparable (normalization and distribution shift)
#######################################################################
raw_heat_hist<-raster("./Input_data/spatial_data/da_n_day_hist.tif")
raw_heat_ssp126<-raster("./Input_data/spatial_data/da_n_day_ssp126.tif")
raw_drought_hist<-raster("./Input_data/spatial_data/SPEI_low10_hist.tif")
raw_drought_ssp126<-raster("./Input_data/spatial_data/SPEI_low10_ssp126.tif")
raw_flood_hist<-raster("./Input_data/spatial_data/SDII_hist.tif")
raw_flood_ssp126<-raster("./Input_data/spatial_data/SDII_ssp126.tif")

raw_heat_hist<-formatRast(raw_heat_hist)
raw_heat_ssp126<-formatRast(raw_heat_ssp126)
raw_drought_hist<-formatRast(raw_drought_hist)
raw_drought_ssp126<-formatRast(raw_drought_ssp126)
raw_flood_hist<-formatRast(raw_flood_hist)
raw_flood_ssp126<-formatRast(raw_flood_ssp126)

##baseline
raw_heat_hist2<-raw_heat_hist
raw_heat_hist2[raw_heat_hist2<1]<-NA
heat_hist_nrm<-(raw_heat_hist2-minValue(raw_heat_hist2))/(maxValue(raw_heat_hist2)-minValue(raw_heat_hist2))
raw_drought_hist2<- -raw_drought_hist
drought_hist_nrm<- (raw_drought_hist2-minValue(raw_drought_hist2))/(maxValue(raw_drought_hist2)-minValue(raw_drought_hist2))
flood_hist_nrm<-(raw_flood_hist-minValue(raw_flood_hist))/(maxValue(raw_flood_hist)-minValue(raw_flood_hist))
##future
raw_heat_ssp126_2<-raw_heat_ssp126
raw_heat_ssp126_2[raw_heat_ssp126_2<1]<-NA
heat_ssp126_nrm<-(raw_heat_ssp126_2-minValue(raw_heat_ssp126_2))/(maxValue(raw_heat_ssp126_2)-minValue(raw_heat_ssp126_2))
raw_drought_ssp126_2 <- -raw_drought_ssp126
drought_ssp126_nrm<- (raw_drought_ssp126_2-minValue(raw_drought_ssp126_2))/(maxValue(raw_drought_ssp126_2)-minValue(raw_drought_ssp126_2))
flood_ssp126_nrm<-(raw_flood_ssp126-minValue(raw_flood_ssp126))/(maxValue(raw_flood_ssp126)-minValue(raw_flood_ssp126))

#plotting the distribution of raw values
tiff('./plots/raw_hazard_distribution_ssp126.tiff', units="in", width=7, height=4, res=500,compression = 'lzw')
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

hist(raw_heat_ssp126_2,xlab='number if days with HI>41°C',ylab='frequency',breaks = "Scott",main='')
title(main = "heat",adj = 0)
abline(v=quantile(raw_heat_ssp126_2, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=30,col='red',lwd=1.5,lty=2)
hist(raw_drought_ssp126_2,xlab='-(SPEI)',ylab='frequency',breaks=80,main='')
title(main = "drought",adj = 0)
abline(v=quantile(raw_drought_ssp126_2, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=1,col='red',lwd=1.5,lty=2)
mtext("Future (SSP126)", side = 3, line = 3, outer = FALSE)
hist(raw_flood_ssp126,xlab='SDII',ylab='frequency',breaks = "Scott",main='')
title(main = "flood",adj = 0)
abline(v=quantile(raw_flood_ssp126, 0.8 , na.rm = TRUE),col='blue',lwd=1.5,lty=1)
abline(v=10,col='red',lwd=1.5,lty=2)

dev.off()

#set the threshold to be considered as 'severe' (class 4), and normalize their values
heat_thresh_nrm_hist<-(30-minValue(raw_heat_hist2))/(maxValue(raw_heat_hist2)-minValue(raw_heat_hist2))
drought_thresh_nrm_hist<- (1-minValue(raw_drought_hist2))/(maxValue(raw_drought_hist2)-minValue(raw_drought_hist2))
flood_thresh_nrm_hist<- (10-minValue(raw_flood_hist))/(maxValue(raw_flood_hist)-minValue(raw_flood_hist))
heat_thresh_nrm_ssp126<-(30-minValue(raw_heat_ssp126))/(maxValue(raw_heat_ssp126)-minValue(raw_heat_ssp126))
drought_thresh_nrm_ssp126<- (1-minValue(raw_drought_ssp126_2))/(maxValue(raw_drought_ssp126_2)-minValue(raw_drought_ssp126_2))
flood_thresh_nrm_ssp126<- (10-minValue(raw_flood_ssp126))/(maxValue(raw_flood_ssp126)-minValue(raw_flood_ssp126))

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
tiff('./plots/hazard_distribution_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
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

df_freq_ssp126<- data.frame(x = values(heat_ssp126_nrm), y = values(drought_ssp126_nrm),z=values(flood_ssp126_nrm)) %>%
  gather(key, value)
tiff('./plots/hazard_distribution_ssp126.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_ssp126, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = heat_thresh_nrm_ssp126, linetype="dashed", 
             color = "red")+
  geom_vline(xintercept = drought_thresh_nrm_ssp126, linetype="dashed", 
             color = "yellow")+
  geom_vline(xintercept = flood_thresh_nrm_ssp126, linetype="dashed", 
             color = "blue")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - future (SSP126)',x='hazard index value',y='density')
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
heat_ssp126_sft<-heat_ssp126_nrm+(0.6-heat_thresh_nrm_ssp126)
drought_ssp126_sft<-drought_ssp126_nrm+(0.6-drought_thresh_nrm_ssp126)
flood_ssp126_sft<-flood_ssp126_nrm+(0.6-flood_thresh_nrm_ssp126)
heat_ssp126_sft[heat_ssp126_sft>=1]<-1
drought_ssp126_sft[drought_ssp126_sft>=1]<-1
flood_ssp126_sft[flood_ssp126_sft>=1]<-1
heat_ssp126_sft[heat_ssp126_sft<=0]<-0
drought_ssp126_sft[drought_ssp126_sft<=0]<-0
flood_ssp126_sft[flood_ssp126_sft<=0]<-0

heat_thresh_sft_ssp126<- heat_thresh_nrm_ssp126 + (0.6-heat_thresh_nrm_ssp126) 
drought_thresh_sft_ssp126<- drought_thresh_nrm_ssp126 + (0.6-drought_thresh_nrm_ssp126) 
flood_thresh_sft_ssp126<- flood_thresh_nrm_ssp126 + (0.6-flood_thresh_nrm_ssp126)

#plot the new shifted distributions
df_freq_hist_sft<- data.frame(x = values(heat_hist_sft), y = values(drought_hist_sft),z=values(flood_hist_sft)) %>%
  gather(key, value)
tiff('./plots/hazard_distribution_shifted_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
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

df_freq_ssp126_sft<- data.frame(x = values(heat_ssp126_sft), y = values(drought_ssp126_sft),z=values(flood_ssp126_sft)) %>%
  gather(key, value)
tiff('./plots/hazard_distribution_shifted_ssp126.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_ssp126_sft, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = heat_thresh_sft_ssp126, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = drought_thresh_sft_ssp126, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = flood_thresh_sft_ssp126, linetype="dashed", 
             color = "black")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - future (SSP126)',x='hazard index value',y='density')
dev.off()

heat_hist_index<-heat_hist_sft
heat_hist_index[raw_heat_hist<1]<-minValue(heat_hist_index)
drought_hist_index<-drought_hist_sft
flood_hist_index<-flood_hist_sft

heat_ssp126_index<-heat_ssp126_sft
heat_ssp126_index[raw_heat_ssp126<1]<-minValue(heat_ssp126_index)
drought_ssp126_index<-drought_ssp126_sft
flood_ssp126_index<-flood_ssp126_sft

heat_change_index<-heat_ssp126_index-heat_hist_index
drought_change_index<-drought_ssp126_index-drought_hist_index
flood_change_index<-flood_ssp126_index-flood_hist_index

#plotting the distribution of the new normalized hazard index
tiff('./plots/hazard_distribution_ssp126.tiff', units="in", width=7, height=4, res=500,compression = 'lzw')
par(mfrow = c(2, 3))
hist(heat_hist_index,xlab='heat index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "heat",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
hist(drought_hist_index,xlab='drought index',ylab='frequency',breaks=75,main='',xlim=c(0,1))
title(main = "drought",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
mtext("Baseline", side = 3, line = 3, outer =FALSE)
hist(flood_hist_index,xlab='flood index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "flood",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
legend('topright',inset=c(-0.1, -0.5), legend=c("class threshold"),
       col='blue', lty=2, cex=0.8,xpd=TRUE)

hist(heat_ssp126_index,xlab='heat index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "heat",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
hist(drought_ssp126_index,xlab='drought index',ylab='frequency',breaks=75,main='',xlim=c(0,1))
title(main = "drought",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
mtext("Future (SSP126)", side = 3, line = 3, outer =FALSE)
hist(flood_ssp126_index,xlab='flood index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "flood",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)


dev.off()


#remove Greenland from the maps
gr<-HoA[which(HoA$gis_name=='Greenland (DNK)'),]
heat_hist_index<-mask(heat_hist_index,gr,inverse=TRUE)
heat_ssp126_index<-mask(heat_ssp126_index,gr,inverse=TRUE)
drought_hist_index<-mask(drought_hist_index,gr,inverse=TRUE)
drought_ssp126_index<-mask(drought_ssp126_index,gr,inverse=TRUE)
flood_hist_index<-mask(flood_hist_index,gr,inverse=TRUE)
flood_ssp126_index<-mask(flood_ssp126_index,gr,inverse=TRUE)
#write the indices into raster files
writeRaster(heat_hist_index, filename = ("./output_data/spatial_data/heat_hist_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(drought_hist_index, filename = ("./output_data/spatial_data/drought_hist_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(flood_hist_index, filename = ("./output_data/spatial_data/flood_hist_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(heat_ssp126_index, filename = ("./output_data/spatial_data/heat_ssp126_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(drought_ssp126_index, filename = ("./output_data/spatial_data/drought_ssp126_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(flood_ssp126_index, filename = ("./output_data/spatial_data/flood_ssp126_nrm.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(heat_change_index, filename = ("./output_data/spatial_data/heat_change_nrm_ssp126.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(drought_change_index, filename = ("./output_data/spatial_data/drought_change_nrm_ssp126.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(flood_change_index, filename = ("./output_data/spatial_data/flood_change_nrm_ssp126.tif"), format = "GTiff",overwrite=TRUE)


##visualizing individual climate hazard indices (current, future, change)
##visualizing heat hazard
heat_hist_tbl <- overlayRast(heat_hist_index)
heat_ssp126_tbl <- overlayRast(heat_ssp126_index)
heat_change_index<-heat_ssp126_index-heat_hist_index
heat_change_tbl<-overlayRast(heat_change_index)
the_min<-min(c(minValue(heat_hist_index),minValue(heat_ssp126_index)))
the_max<-max(c(maxValue(heat_hist_index),maxValue(heat_ssp126_index)))
col_map<-rev(sequential_hcl(27, palette = "Reds3"))
col_map<-col_map[2:length(col_map)]
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster("./plots/heat_hist_index.tiff",color_bar_min,color_bar_max,"heat hazard - baseline",'heat index',col_map,the_ticks,the_labels,heat_hist_tbl)
Visu_raster("./plots/heat_ssp126_index.tiff",color_bar_min,color_bar_max,"heat hazard - future (SSP126)",'heat index',col_map,the_ticks,the_labels,heat_ssp126_tbl)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster("./plots/heat_change_index_ssp126.tiff",minValue(heat_change_index),maxValue(heat_change_index),"heat hazard - change (SSP126)",'heat change',col_map,waiver(),waiver(),heat_change_tbl)

##visualizing drought hazard
drought_hist_tbl <- overlayRast(drought_hist_index)
drought_ssp126_tbl <- overlayRast(drought_ssp126_index)
drought_change_index<-drought_ssp126_index-drought_hist_index
drought_change_tbl<-overlayRast(drought_change_index)
the_min<-min(c(minValue(drought_hist_index),minValue(drought_ssp126_index)))
the_max<-max(c(maxValue(drought_hist_index),maxValue(drought_ssp126_index)))
col_map<-rev(sequential_hcl(40, palette = "YlOrBr"))
col_map<-c(rep(col_map[1],3),rep(col_map[2],4),rep(col_map[3],5),col_map[4:length(col_map)])
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster("./plots/drought_hist_index.tiff",color_bar_min,color_bar_max,"drought hazard - baseline",'drought index',col_map,the_ticks,the_labels,drought_hist_tbl)
Visu_raster("./plots/drought_ssp126_index.tiff",color_bar_min,color_bar_max,"drought hazard - future (SSP126)",'drought index',col_map,the_ticks,the_labels,drought_ssp126_tbl)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster("./plots/drought_change_index_ssp126.tiff",minValue(drought_change_index),maxValue(drought_change_index),"drought hazard - change (SSP126)",'drought change',col_map,waiver(),waiver(),drought_change_tbl)

##visualizing flood hazard
flood_hist_tbl <- overlayRast(flood_hist_index)
flood_ssp126_tbl <- overlayRast(flood_ssp126_index)
flood_change_index<-flood_ssp126_index-flood_hist_index
flood_change_tbl<-overlayRast(flood_change_index)
the_min<-min(c(minValue(flood_hist_index),minValue(flood_ssp126_index)))
the_max<-max(c(maxValue(flood_hist_index),maxValue(flood_ssp126_index)))
col_map<-rev(sequential_hcl(40, palette = "Blues3"))
col_map<-c(rep(col_map[1],3),rep(col_map[2],4),rep(col_map[3],5),col_map[4:length(col_map)])
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster("./plots/flood_hist_index.tiff",color_bar_min,color_bar_max,"flood hazard - baseline",'flood index',col_map,the_ticks,the_labels,flood_hist_tbl)
Visu_raster("./plots/flood_ssp126_index.tiff",color_bar_min,color_bar_max,"flood hazard - future (SSP126)",'flood index',col_map,the_ticks,the_labels,flood_ssp126_tbl)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster("./plots/flood_change_index_ssp126.tiff",minValue(flood_change_index),maxValue(flood_change_index),"flood hazard - change (SSP126)",'flood change',col_map,waiver(),waiver(),flood_change_tbl)
