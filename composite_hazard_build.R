# This script has two parts (1) building the composite hazard index based on the individual hazards indices
# (2) investigate on the individual hazard that is mostly driving the high compounding effects and how these hazards are distributed over regions 

packages <- c("viridis", "pacman", "DescTools", "raster", "rgdal", "tidyverse", "tibble", "ggplot2", "sf", "scales",
              "RColorBrewer","cartography", "tidyr", "openxlsx")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

HoA <- st_read("./input_data/spatial_data/world_shape2.geojson")
##############################################
#defining all functions and plot formatting
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

##reading individual indices from files
heat_hist_index<-raster(("./output_data/spatial_data/heat_hist_nrm.tif"))
heat_ssp_index<-raster(("./output_data/spatial_data/heat_ssp_nrm.tif"))
drought_hist_index<-raster(("./output_data/spatial_data/drought_hist_nrm.tif"))
drought_ssp_index<-raster(("./output_data/spatial_data/drought_ssp_nrm.tif"))
flood_hist_index<-raster(("./output_data/spatial_data/flood_hist_nrm.tif"))
flood_ssp_index<-raster(("./output_data/spatial_data/flood_ssp_nrm.tif"))


#########################part 1: building the composite hazard index ############
#################################################################################
#building composite index
composite_hist_index<- -(1/(sqrt((heat_hist_index^2)+(drought_hist_index^2)+(flood_hist_index^2))))
composite_ssp_index<- -(1/(sqrt((heat_ssp_index^2)+(drought_ssp_index^2)+(flood_ssp_index^2))))

mi<-min(c(minValue(composite_hist_index),minValue(composite_ssp_index)))
ma<-max(c(maxValue(composite_hist_index),maxValue(composite_ssp_index)))
composite_hist_index<-(composite_hist_index-mi)/(ma-mi)
composite_ssp_index<-(composite_ssp_index-mi)/(ma-mi)

composite_change_index<-composite_ssp_index-composite_hist_index

#remove Greenland from the maps
gr<-HoA[which(HoA$gis_name=='Greenland (DNK)'),]
composite_hist_index<-mask(composite_hist_index,gr,inverse=TRUE)
composite_ssp_index<-mask(composite_ssp_index,gr,inverse=TRUE)
composite_change_index<-mask(composite_change_index,gr,inverse=TRUE)

#plotting distribution of compound hazard
composite_hist_tbl<- overlayRast(composite_hist_index)
composite_ssp_tbl <- overlayRast(composite_ssp_index)
composite_change_tbl <- overlayRast(composite_change_index)

the_min<-min(c(minValue(composite_hist_index),minValue(composite_ssp_index)))
the_max<-max(c(maxValue(composite_hist_index),maxValue(composite_ssp_index)))
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c(the_ticks[1],rep('',length(the_ticks)-2),the_ticks[length(the_ticks)])
col_map<-rev(brewer.pal(11,'Spectral'))
Visu_raster("./plots/composite_hist.tiff",color_bar_min,color_bar_max,"composite hazard - baseline",'composite index',col_map,the_ticks,the_labels,composite_hist_tbl)
Visu_raster("./plots/composite_ssp.tiff",color_bar_min,color_bar_max,"composite hazard - future",'composite index',col_map,the_ticks,the_labels,composite_ssp_tbl)
Visu_raster("./plots/composite_change.tiff",minValue(composite_change_index),maxValue(composite_change_index),"composite hazard - change",'composite change',col_map,waiver(),waiver(),composite_change_tbl)

#saving the compound climate hazard maps into rasters
writeRaster(composite_hist_index, filename = ("./output_data/spatial_data/composite_hist_index.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(composite_ssp_index, filename = ("./output_data/spatial_data/composite_ssp_index.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(composite_change_index, filename = ("./output_data/spatial_data/composite_change_index.tif"), format = "GTiff",overwrite=TRUE)

# build the 5 different severity categories and plot the resulting maps
brk  <- seq(0,1,0.2)
composite_hist_class <- raster::cut(composite_hist_index, breaks=brk) 
composite_hist_class_tbl<-overlayRast(composite_hist_class)

composite_ssp_class <- raster::cut(composite_ssp_index, breaks=brk) 
composite_ssp_class_tbl<-overlayRast(composite_ssp_class)

spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
Visu_raster_discrete("./plots/composite_hist_class.tiff","composite hazard - baseline",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_hist_class_tbl)
Visu_raster_discrete("./plots/composite_ssp_class.tiff","composite hazard - future",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_ssp_class_tbl)

#saving the class maps into rasters for future use
writeRaster(composite_hist_class, filename = ("./output_data/spatial_data/composite_hist_class.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(composite_ssp_class, filename = ("./output_data/spatial_data/composite_ssp_class.tif"), format = "GTiff",overwrite=TRUE)

#average the composite hazard per country 
sp_data_hist<-raster::extract(composite_hist_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_hist<-composite_hist_index
values(empty_rast_hist)<-NA
the_mean_hist<-sp_data_hist[dim(sp_data_hist)[2]]
hist_avrg<-rasterize(the_mean_hist, empty_rast_hist, field = the_mean_hist@data[,1], update = TRUE)
hist_avrg_tbl<-overlayRast(hist_avrg)

sp_data_ssp<-raster::extract(composite_ssp_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_ssp<-composite_ssp_index
values(empty_rast_ssp)<-NA
the_mean_ssp<-sp_data_ssp[dim(sp_data_ssp)[2]]
ssp_avrg<-rasterize(the_mean_ssp, empty_rast_ssp, field = the_mean_ssp@data[,1], update = TRUE)
ssp_avrg_tbl<-overlayRast(ssp_avrg)

col_map<-rev(brewer.pal(11,'Spectral'))
Visu_raster("./plots/composite_hist_avrg.tiff",0,1,"composite hazard - baseline",'composite index',col_map,waiver(),waiver(),hist_avrg_tbl)
Visu_raster("./plots/composite_ssp_avrg.tiff",0,1,"composite hazard - future",'composite index',col_map,waiver(),waiver(),ssp_avrg_tbl)

#classify each country within the 5 categories
brk  <- seq(0,1,0.2)
hist_class_avrg <- raster::cut(hist_avrg, breaks=brk)
hist_class_avrg_tbl<-overlayRast(hist_class_avrg)

ssp_class_avrg <- raster::cut(ssp_avrg, breaks=brk)
ssp_class_avrg_tbl<-overlayRast(ssp_class_avrg)


spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
Visu_raster_discrete("./plots/composite_hist_class_avrg.tiff","composite hazard - baseline",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),hist_class_avrg_tbl)
Visu_raster_discrete("./plots/composite_ssp_class_avrg.tiff","composite hazard - future",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),ssp_class_avrg_tbl)

#saving the class per country maps into rasters for future use
writeRaster(hist_class_avrg, filename = ("./output_data/spatial_data/composite_hist_class_avrg.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(ssp_class_avrg, filename = ("./output_data/spatial_data/composite_ssp_class_avrg.tif"), format = "GTiff",overwrite=TRUE)


#######################  part 2: investigate on the individual hazard that is mostly driving the high compounding effects 
#          and how these hazards are distributed over regions 
######################

#plot the individual hazards under the high, severe and extreme compound hazard 

comp_hist_from_high<-composite_hist_class
comp_hist_from_high[comp_hist_from_high<=2]<-NA
comp_ssp_from_high<-composite_ssp_class
comp_ssp_from_high[comp_ssp_from_high<=2]<-NA

the_point<-0.4

high_ind_hist <- comp_hist_from_high
high_ind_hist[,]<-0
high_ind_hist[heat_hist_index>=the_point & drought_hist_index<the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)]<-1
high_ind_hist[heat_hist_index<the_point & drought_hist_index>=the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)]<-2
high_ind_hist[heat_hist_index<the_point & drought_hist_index<the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-3
high_ind_hist[heat_hist_index>=the_point & drought_hist_index>=the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)]<-4
high_ind_hist[heat_hist_index>=the_point & drought_hist_index<the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-5
high_ind_hist[heat_hist_index<the_point & drought_hist_index>=the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-6
high_ind_hist[heat_hist_index>=the_point & drought_hist_index>=the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-7
high_ind_hist[heat_hist_index<the_point & drought_hist_index<the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)  ]<-7

high_ind_hist_tbl<-overlayRast(high_ind_hist)


high_ind_ssp <- comp_ssp_from_high
high_ind_ssp[,]<-0
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index<the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)]<-1
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index>=the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)]<-2
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index<the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-3
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index>=the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)]<-4
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index<the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-5
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index>=the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-6
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index>=the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-7
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index<the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)  ]<-7

high_ind_ssp_tbl<-overlayRast(high_ind_ssp)

my_pals=c('0'='white','1'='orangered1','2'='goldenrod3','3'='steelblue3','4'='saddlebrown','5'='olivedrab1','6'='purple','7'='grey')

Visu_raster_discrete("./plots/high_to_severe_compound_hist.tiff","high to severe hazard - baseline",'hazard type',my_pals,c('0','1','2','3','4','5','6','7'),c('0','1','2','3','4','5','6','7'),c('','heat','drought','flood','heat and drougth','heat and flood','drought and flood','all'),high_ind_hist_tbl)
Visu_raster_discrete("./plots/high_to_severe_compound_ssp.tiff","high to severe hazard - future",'hazard type',my_pals,c('0','1','2','3','4','5','6','7'),c('0','1','2','3','4','5','6','7'),c('','heat','drought','flood','heat and drougth','heat and flood','drought and flood','all'),high_ind_ssp_tbl)

#plot the distribution of underlaying driving hazard per region
##baseline
regional_files<-c('AMERICAs_15m_UNHCR.gpkg',"ASIA_15m_UNHCR.gpkg","East_and_horn_of_Africa_shapefile.gpkg",
                  "Europe_new.gpkg","MENA_polygon_2.gpkg","Southern_Africa_Ash.gpkg","West_and_central_new.gpkg")

all_region<-c('Americas','Asia','East and\nHorn of Africa','Europe','MENA','Southern\nAfrica','West and\nCentral Africa')
all_types=c('heat','drought','flood','heat and drougth','heat and flood','drought and flood','all')

percentage=c()
type=c()
region=c()
for (count_region in 1:length(regional_files)){
  polygon_rst <- st_read(paste0("./input_data/spatial_data/",regional_files[count_region]))
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
    
    percentage<-c(percentage,proportion_type*100)
    type<-c(type,all_types[type_ind])
    region<-c(region,all_region[count_region])
  }
}
all_colors=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3','heat and drougth'='saddlebrown','heat and flood'='olivedrab1','drought and flood'='purple','all'='grey')
the_df<-data.frame(percentage,type,region)
the_df$region<-factor(the_df$region)
the_df2<-the_df[the_df$percentage!=0,]
tiff('./plots/from_high_per_region_hist.tiff', units="in", width=11, height=6, res=500,compression = 'lzw')
ggplot(the_df,aes(x=region,y=percentage,fill=type))+
  geom_bar(stat='identity')+
  labs(y='area under high, severe and extreme (%)')+
  scale_fill_manual(values=all_colors,breaks=unique(the_df2$type),labels=all_types[match(unique(the_df2$type),all_types)])+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  ylim(0,100)+
  ggtitle('baseline')+
  theme(axis.text=element_text(size=16),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.1,size=25),legend.text = element_text(size = 15),legend.title = element_text(size = 18))
dev.off()

#future
percentage=c()
type=c()
region=c()
for (count_region in 1:length(regional_files)){
  polygon_rst <- st_read(paste0("./input_data/spatial_data/",regional_files[count_region]))
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
    
    percentage<-c(percentage,proportion_type*100)
    type<-c(type,all_types[type_ind])
    region<-c(region,all_region[count_region])
  }
}
the_df<-data.frame(percentage,type,region)
the_df$region<-factor(the_df$region)
the_df2<-the_df[the_df$percentage!=0,]

tiff('./plots/from_high_per_region_ssp.tiff', units="in", width=11, height=6, res=500,compression = 'lzw')
ggplot(the_df,aes(x=region,y=percentage,fill=type))+
  geom_bar(stat='identity')+
  labs(y='area under high, severe and extreme (%)')+
  scale_fill_manual(values=all_colors,breaks=unique(the_df2$type),labels=all_types[match(unique(the_df2$type),all_types)])+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  ylim(0,100)+
  ggtitle('future')+
  theme(axis.text=element_text(size=16),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.1,size=25),legend.text = element_text(size = 15),legend.title = element_text(size = 18))

dev.off()

