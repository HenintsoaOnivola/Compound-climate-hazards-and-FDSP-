#This script looks at the exposure of FDSP population to each individual hazard and their specific class

packages <- c("viridis", "pacman", "DescTools", "raster", "rgdal", "tidyverse", "tibble", "ggplot2", "sf", "scales",
              "RColorBrewer", "tidyr", "spatstat",'maptools','readxl','readr','reshape2','cowplot')

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

HoA <- st_read("./input_data/spatial_data/world_shape2.geojson")


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

#read hazard class indices
###baseline
heat_hist_index<-raster("./output_data/spatial_data/heat_hist_nrm.tif")
heat_data_hist<-raster::extract(heat_hist_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_hist<-heat_hist_index
values(empty_rast_hist)<-NA
the_mean_hist<-heat_data_hist[dim(heat_data_hist)[2]]
heat_hist_avrg<-rasterize(the_mean_hist, empty_rast_hist, field = the_mean_hist@data[,1], update = TRUE)
heat_hist_class_avrg <- raster::cut(heat_hist_avrg, breaks=seq(0,1,0.2))

drought_hist_index<-raster("./output_data/spatial_data/drought_hist_nrm.tif")
drought_hist_class <- raster::cut(drought_hist_index, breaks=seq(0,1,0.2))
drought_data_hist<-raster::extract(drought_hist_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_hist<-drought_hist_index
values(empty_rast_hist)<-NA
the_mean_hist<-drought_data_hist[dim(drought_data_hist)[2]]
drought_hist_avrg<-rasterize(the_mean_hist, empty_rast_hist, field = the_mean_hist@data[,1], update = TRUE)
drought_hist_class_avrg <- raster::cut(drought_hist_avrg, breaks=seq(0,1,0.2))

flood_hist_index<-raster("./output_data/spatial_data/flood_hist_nrm.tif")
flood_hist_class <- raster::cut(flood_hist_index, breaks=seq(0,1,0.2))
flood_data_hist<-raster::extract(flood_hist_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_hist<-flood_hist_index
values(empty_rast_hist)<-NA
the_mean_hist<-flood_data_hist[dim(flood_data_hist)[2]]
flood_hist_avrg<-rasterize(the_mean_hist, empty_rast_hist, field = the_mean_hist@data[,1], update = TRUE)
flood_hist_class_avrg <- raster::cut(flood_hist_avrg, breaks=seq(0,1,0.2))

####future
heat_ssp_index<-raster("./output_data/spatial_data/heat_ssp_nrm.tif")
heat_data_ssp<-raster::extract(heat_ssp_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_ssp<-heat_ssp_index
values(empty_rast_ssp)<-NA
the_mean_ssp<-heat_data_ssp[dim(heat_data_ssp)[2]]
heat_ssp_avrg<-rasterize(the_mean_ssp, empty_rast_ssp, field = the_mean_ssp@data[,1], update = TRUE)
heat_ssp_class_avrg <- raster::cut(heat_ssp_avrg, breaks=seq(0,1,0.2))

drought_ssp_index<-raster("./output_data/spatial_data/drought_ssp_nrm.tif")
drought_ssp_class <- raster::cut(drought_ssp_index, breaks=seq(0,1,0.2))
drought_data_ssp<-raster::extract(drought_ssp_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_ssp<-drought_ssp_index
values(empty_rast_ssp)<-NA
the_mean_ssp<-drought_data_ssp[dim(drought_data_ssp)[2]]
drought_ssp_avrg<-rasterize(the_mean_ssp, empty_rast_ssp, field = the_mean_ssp@data[,1], update = TRUE)
drought_ssp_class_avrg <- raster::cut(drought_ssp_avrg, breaks=seq(0,1,0.2))

flood_ssp_index<-raster("./output_data/spatial_data/flood_ssp_nrm.tif")
flood_ssp_class <- raster::cut(flood_ssp_index, breaks=seq(0,1,0.2))
flood_data_ssp<-raster::extract(flood_ssp_index, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_ssp<-flood_ssp_index
values(empty_rast_ssp)<-NA
the_mean_ssp<-flood_data_ssp[dim(flood_data_ssp)[2]]
flood_ssp_avrg<-rasterize(the_mean_ssp, empty_rast_ssp, field = the_mean_ssp@data[,1], update = TRUE)
flood_ssp_class_avrg <- raster::cut(flood_ssp_avrg, breaks=seq(0,1,0.2))

#read FDSP data
PoC<-read_csv("./input_data/FDSP_population.csv",skip=14)
our_PoC<-PoC[PoC$Year==2022,]
PoC_data<-data.frame(our_PoC)
PoC_data<-cbind(country_name=PoC_data[,4],country_code=PoC_data[,5],PoC_data[,c(6:10,12)])
PoC_data[PoC_data=='-']<-NA
PoC_data[,6]<-as.numeric(PoC_data[,6])
PoC_number<-rowSums(PoC_data[,3:8],na.rm=TRUE)
PoC_data<-data.frame(country_name=PoC_data[,1],country_code=PoC_data[,2],PoC_number)
PoC_data<-aggregate(PoC_number~(country_code+country_name),data=PoC_data,FUN=sum)

#add hazard class for each country within the FDSP data 
#########################

hist_class_avrg<-raster("./output_data/spatial_data/composite_hist_class_avrg.tif")
ssp_class_avrg<-raster("./output_data/spatial_data/composite_ssp_class_avrg.tif")

comp_class_hist<-c()
comp_class_ssp<-c()
heat_class_hist<-c()
heat_class_ssp<-c()
drought_class_hist<-c()
drought_class_ssp<-c()
flood_class_hist<-c()
flood_class_ssp<-c()

for (i in 1:nrow(PoC_data)){
  the_ind<-which(HoA$iso3==(PoC_data$country_code[i]))
  if (length(the_ind)==0){
    the_ind<-grep((PoC_data$country_name)[i],HoA$gis_name)
  }
  if (length(the_ind)!=0){
    the_ind1<-the_ind[1]
    the_country_code<-(HoA$iso3)[the_ind1]
    
    #hazard class per country
    all_haz_hist<-unlist(raster::extract(hist_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz1<- (as.integer(names(sort(table(all_haz_hist),decreasing=TRUE)[1])))
    heat_haz_hist<-unlist(raster::extract(heat_hist_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz2<- (as.integer(names(sort(table(heat_haz_hist),decreasing=TRUE)[1])))
    drought_haz_hist<-unlist(raster::extract(drought_hist_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz3<- (as.integer(names(sort(table(drought_haz_hist),decreasing=TRUE)[1])))
    flood_haz_hist<-unlist(raster::extract(flood_hist_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz4<- (as.integer(names(sort(table(flood_haz_hist),decreasing=TRUE)[1])))
    if (length(the_haz1)==0){
      comp_class_hist<-c(comp_class_hist,NA)
      heat_class_hist<-c(heat_class_hist,NA)
      drought_class_hist<-c(drought_class_hist,NA)
      flood_class_hist<-c(flood_class_hist,NA)
      
    }else{
      comp_class_hist<-c(comp_class_hist,the_haz1)
      heat_class_hist<-c(heat_class_hist,the_haz2)
      drought_class_hist<-c(drought_class_hist,the_haz3)
      flood_class_hist<-c(flood_class_hist,the_haz4)
    }
    
    all_haz_ssp<-unlist(raster::extract(ssp_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz1<- (as.integer(names(sort(table(all_haz_ssp),decreasing=TRUE)[1])))
    heat_haz_ssp<-unlist(raster::extract(heat_ssp_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz2<- (as.integer(names(sort(table(heat_haz_ssp),decreasing=TRUE)[1])))
    drought_haz_ssp<-unlist(raster::extract(drought_ssp_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz3<- (as.integer(names(sort(table(drought_haz_ssp),decreasing=TRUE)[1])))
    flood_haz_ssp<-unlist(raster::extract(flood_ssp_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz4<- (as.integer(names(sort(table(flood_haz_ssp),decreasing=TRUE)[1])))
    if (length(the_haz1)==0){
      comp_class_ssp<-c(comp_class_ssp,NA)
      heat_class_ssp<-c(heat_class_ssp,NA)
      drought_class_ssp<-c(drought_class_ssp,NA)
      flood_class_ssp<-c(flood_class_ssp,NA)
      
    }else{
      comp_class_ssp<-c(comp_class_ssp,the_haz1)
      heat_class_ssp<-c(heat_class_ssp,the_haz2)
      drought_class_ssp<-c(drought_class_ssp,the_haz3)
      flood_class_ssp<-c(flood_class_ssp,the_haz4)
    }

  }else{
    #look in the other dictionary for region
    comp_class_hist<-c(comp_class_hist,NA)
    heat_class_hist<-c(heat_class_hist,NA)
    drought_class_hist<-c(drought_class_hist,NA)
    flood_class_hist<-c(flood_class_hist,NA)
    comp_class_ssp<-c(comp_class_ssp,NA)
    heat_class_ssp<-c(heat_class_ssp,NA)
    drought_class_ssp<-c(drought_class_ssp,NA)
    flood_class_ssp<-c(flood_class_ssp,NA)
  }
}
PoC_data<-cbind(PoC_data,comp_class_hist,comp_class_ssp,heat_class_hist,heat_class_ssp,drought_class_hist,drought_class_ssp,flood_class_hist,flood_class_ssp)
saveRDS(PoC_data,'./output_data/FDSP_population_data_frame_ind_hazard.Rds')
PoC_data<-readRDS('./output_data/FDSP_population_data_frame_ind_hazard.Rds')

#Plot FDSP numbers per hazard class
exposure_hist<-data.frame(hazard.class=rep(1:5,4),PoC.number=NA,hazard=c(rep('heat',5),rep('drought',5),rep('flood',5),rep('compound',5)))
exposure_ssp<-data.frame(hazard.class=rep(1:5,4),PoC.number=NA,hazard=c(rep('heat',5),rep('drought',5),rep('flood',5),rep('compound',5)))
for (h in c('heat','drought','flood','compound')){
  for (cl in 1:5){
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='heat',2]<-PoC_data$PoC_number[which((PoC_data$heat_class_hist)==cl)]%>%sum
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='drought',2]<-PoC_data$PoC_number[which((PoC_data$drought_class_hist)==cl)]%>%sum
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='flood',2]<-PoC_data$PoC_number[which((PoC_data$flood_class_hist)==cl)]%>%sum
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='compound',2]<-PoC_data$PoC_number[which((PoC_data$comp_class_hist)==cl)]%>%sum
    
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='heat',2]<-PoC_data$PoC_number[which((PoC_data$heat_class_ssp)==cl)]%>%sum
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='drought',2]<-PoC_data$PoC_number[which((PoC_data$drought_class_ssp)==cl)]%>%sum
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='flood',2]<-PoC_data$PoC_number[which((PoC_data$flood_class_ssp)==cl)]%>%sum
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='compound',2]<-PoC_data$PoC_number[which((PoC_data$comp_class_ssp)==cl)]%>%sum
  }
}

spect_pal<-rev(brewer.pal(5,'Spectral'))
all_colors=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
all_class=names(all_colors)

exposure_hist$hazard<-factor(exposure_hist$hazard)
exposure_hist$hazard.class<-factor(exposure_hist$hazard.class)
the_df<-exposure_hist[exposure_hist$PoC.number!=0,]
the_df$PoC.number<-the_df$PoC.number/1e6

tiff('./plots/with_FDSP/exposure_per_class.tiff', units="in", width=8, height=4, res=500,compression = 'lzw')
ggplot(the_df,aes(x=hazard,y=PoC.number,fill=hazard.class))+
  geom_bar(stat='identity',color = "black",size=0.2,position = position_stack(reverse = TRUE))+
  labs(y='FDSP number (millions)')+
  scale_fill_manual(values=all_colors,name='hazard class')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  theme(axis.text=element_text(size=13),axis.title = element_text(size = 15),legend.text = element_text(size = 13),legend.title = element_text(size = 15))+
  coord_flip()
dev.off()

exposure_ssp$hazard<-factor(exposure_ssp$hazard)
exposure_ssp$hazard.class<-factor(exposure_ssp$hazard.class)
the_df2<-exposure_ssp[exposure_ssp$PoC.number!=0,]
the_df2$PoC.number<-the_df2$PoC.number/1e6

tiff('./plots/with_FDSP/exposure_per_class_ssp.tiff', units="in", width=8, height=4, res=500,compression = 'lzw')
ggplot(the_df2,aes(x=hazard,y=PoC.number,fill=hazard.class))+
  geom_bar(stat='identity',color = "black",size=0.2,position = position_stack(reverse = TRUE))+
  labs(y='FDSP number (millions)')+
  scale_fill_manual(values=all_colors,name='hazard class')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  theme(axis.text=element_text(size=13),axis.title = element_text(size = 15),legend.text = element_text(size = 13),legend.title = element_text(size = 15))+
  coord_flip()
dev.off()

#plot bars wrt class
hist_df<-exposure_hist[exposure_hist$hazard!='compound',]
hist_df$PoC.number<-hist_df$PoC.number/1e6
hist_df$hazard.class<-factor(hist_df$hazard.class)
hist_df$hazard<-factor(hist_df$hazard,levels=c('heat','drought','flood'))

ssp_df<-exposure_ssp[exposure_ssp$hazard!='compound',]
ssp_df$PoC.number<-ssp_df$PoC.number/1e6
ssp_df$hazard.class<-factor(ssp_df$hazard.class)
ssp_df$hazard<-factor(ssp_df$hazard,levels=c('heat','drought','flood'))

my_pals=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3')

p1<-ggplot() + geom_bar(data = hist_df, aes(x = hazard.class, y = PoC.number, fill = hazard), position = "dodge", stat = "identity",width=0.8)+
  labs(x='hazard class', y='FSDSP number (million)')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 12),legend.text = element_text(size = 9),legend.title = element_text(size = 12),legend.key.size = unit(.5, "cm"),
        legend.position = c(0.9, 0.8),axis.text.y = element_text(angle=45,vjust = 0.5, hjust=0.5))+
  scale_x_discrete(labels=c('low','moderate','high','severe','extreme'),breaks=c("1","2","3",'4','5'))+
  scale_fill_manual( values = my_pals)

p1<-p1+coord_flip()
tiff('./plots/with_FDSP/bar_exposure_per_class_hist.tiff',units="in", width=7, height=4, res=500,compression = 'lzw')
p1
dev.off()   


p2<-ggplot() + geom_bar(data = ssp_df, aes(x = hazard.class, y = PoC.number, fill = hazard), position = "dodge", stat = "identity",width=0.8)+
  labs(x='hazard class', y='FDSP number (million)')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 12),legend.text = element_text(size = 9),legend.title = element_text(size = 12),legend.key.size = unit(.5, "cm"),
        legend.position = c(0.9, 0.8),axis.text.y = element_text(angle=45,vjust = 0.5, hjust=0.5))+
  scale_x_discrete(labels=c('low','moderate','high','severe','extreme'),breaks=c("1","2","3",'4','5'))+
  scale_fill_manual( values = my_pals)
p2<-p2+coord_flip()
tiff('./plots/with_FDSP/bar_exposure_per_class_ssp.tiff',units="in", width=7, height=4, res=500,compression = 'lzw')
p2
dev.off() 

tiff('./plots/with_FDSP/bar_exposure_per_class2.tiff',units="in", width=7, height=4, res=500,compression = 'lzw')
prow<-plot_grid(p1+theme(legend.position="none"),p2+theme(legend.position="none"),labels=c('A','B'),align = 'vh',hjust = -1,nrow = 1)
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow, legend, rel_widths = c(2.7, .4))
dev.off()








