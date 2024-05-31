#This script contain the statistical analyses with respect to the exposure of FDSP to each class of climate hazards
#with and without governance taken into account

packages <- c("viridis", "pacman", "DescTools", "raster", "rgdal", "tidyverse", "tibble", "ggplot2", "sf", "scales",
              "RColorBrewer","cartography", "tidyr", "ggrepel",'spatstat','maptools','readxl','readr','ggnewscale')

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
        legend.title = element_text(size=13,face='bold'),
        legend.position = 'bottom',
        legend.box="vertical", 
        legend.margin=margin(),
        legend.key.size = unit(1.5,"line"),
        legend.text=element_text(size=12)
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

Visu_raster_discrete_with_points<-function (output_name,title,bar_title,color_palette,the_lim,all_ticks,all_labels,the_tbl,the_points){
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  
  print(
    ggplot() +
      geom_tile(data = the_tbl, aes(x = x, y = y, fill = factor(value))) +
      geom_sf(data = HoA, fill = NA,lwd=0.5) +
      scale_fill_manual(values=color_palette, name=bar_title,limits=the_lim,breaks=all_ticks,labels=all_labels) +
      new_scale_fill()+  
      geom_sf(data=the_points,size=0.5,aes(fill=new_col))+
      scale_fill_manual(name='',values=c('black'),labels='Location of forcibly displaced and stateless people')+
      guides(fill = guide_legend(override.aes = list(size = 2))) +
      labs(title = title) +
      theme_MAP)
  dev.off()
}

scale_governance<-function(the_value,BTI_min,BTI_max,scaled_min,scaled_max)
{
  (((the_value-BTI_min)/(BTI_min-BTI_max))*(scaled_min-scaled_max))+scaled_min
}

#read the previously built composite hazard indices from rasters
hist_avrg<-raster("./output_data/spatial_data/composite_hist_avrg.tif")
ssp_avrg<-raster("./output_data/spatial_data/composite_ssp_avrg.tif")
composite_hist_index<-raster("./output_data/spatial_data/composite_hist_index.tif")
composite_ssp_index<-raster("./output_data/spatial_data/composite_ssp_index.tif")

composite_hist_class<-raster("./output_data/spatial_data/composite_hist_class.tif")
composite_ssp_class<-raster("./output_data/spatial_data/composite_ssp_class.tif")

brk  <- seq(0,1,0.2)
hist_class_avrg <- raster::cut(hist_avrg, breaks=brk)
ssp_class_avrg <- raster::cut(ssp_avrg, breaks=brk)

#read the governance index per country
BTI_2022<-read_excel("./input_data/BTI_2006-2022_Scores.xlsx",sheet='BTI 2022')
country<-as.data.frame(BTI_2022[1:137,1])
names(country)<-'country'
BTI_ind_2022<-as.data.frame(BTI_2022[,54])
names(BTI_ind_2022)<-'GI_2022'
BTI_2020<-read_excel("./input_data/BTI_2006-2022_Scores.xlsx",sheet='BTI 2020')
BTI_ind_2020<-as.data.frame(BTI_2020[,54])
names(BTI_ind_2020)<-'GI_2020'
BTI_2018<-read_excel("./input_data/BTI_2006-2022_Scores.xlsx",sheet='BTI 2018')
BTI_ind_2018<-as.data.frame(BTI_2018[,54])
names(BTI_ind_2018)<-'GI_2018'

BTI<-cbind(BTI_ind_2022,BTI_ind_2020,BTI_ind_2018)
BTI[BTI=='-']<-NA
BTI<-apply(BTI, 2, FUN=as.numeric)
BTI<-apply(BTI,1,FUN=mean,na.rm = TRUE)
gov_data<-data.frame(country,BTI)

#scale the governance index to account for mis-governance and accommodate the categories defined from BTI index
neg_gov=-gov_data$BTI
scale_frame<-matrix(c(-10,-7,-5.6,-4.3,-3,-7,-5.6,-4.3,-3,-1, 0,0.2,0.4,0.6,0.8,0.2,0.4,0.6,0.8,1),nrow=5,ncol=4)
mis_gov<-rep(NA,length(neg_gov))
for (i in neg_gov){
  ii=which((i<=scale_frame[,2]) & (i>scale_frame[,1]))
  mis_gov[neg_gov==i]<-scale_governance(i,scale_frame[ii,1],scale_frame[ii,2],scale_frame[ii,3],scale_frame[ii,4])
}
gov_data<-cbind(gov_data,misgovernance=mis_gov)

gov_ind<-NA
misgov_ind<-NA
HoA<-cbind(HoA,gov_ind,misgov_ind)
alt_country_name<-c()
c_code<-c()
for (each_c in (gov_data$country)){
  the_ind<-grep(each_c,HoA$gis_name)
  if (length(the_ind)==0){
    var=FALSE
    if (each_c=='Turkey'){
      the_ind<-grep('Türkiye',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Türkiye')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Vietnam'){
      the_ind<-grep('Viet Nam',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Viet Nam')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='South Korea'){
      the_ind<-grep('Republic of Korea',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Republic of Korea')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='North Korea'){
      the_ind<-grep('Democratic People\'s Rep. of Korea',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Democratic People\'s Rep. of Korea')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Congo, DR'){
      the_ind<-grep('Democratic Republic of the Congo',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Democratic Republic of the Congo')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Congo, Rep.'){
      the_ind<-grep('Republic of the Congo',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Republic of the Congo')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Laos'){
      the_ind<-grep('Lao People\'s Democratic Republic',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Lao People\'s Democratic Republic')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Czech Republic'){
      the_ind<-grep('Czechia',HoA$gis_name)
      alt_country_name<-c(alt_country_name,'Czechia')
      c_code<-c(c_code,(HoA$iso3[the_ind])[1])
      var=TRUE
    }
    if (var==FALSE){
      alt_country_name<-c(alt_country_name,each_c)
      c_code<-c(c_code,NA)
    }

  }
  else{
    alt_country_name<-c(alt_country_name,HoA$gis_name[the_ind[1]])
    c_code<-c(c_code,(HoA$iso3[the_ind])[1])
  }
    
  HoA$gov_ind[the_ind]<-gov_data$BTI[gov_data$country==each_c]
  HoA$misgov_ind[the_ind]<-gov_data$misgovernance[gov_data$country==each_c]
}
gov_data<-cbind(gov_data,alt_country_name,c_code)
gov_data<-data.frame(gov_data$country,gov_data$alt_country_name,gov_data$c_code, gov_data$BTI,gov_data$misgovernance)
names(gov_data)<-c('country','country_name_UNHCR','country_code','governance','misgovernance')

rast_base<-hist_avrg
values(rast_base)<-NA
Bgov_rast<-rasterize(HoA,rast_base,field='misgov_ind')

#build the new risk index of hazard and misgovernance
risk_hist<-(composite_hist_index+Bgov_rast)/2
risk_ssp<-(composite_ssp_index+Bgov_rast)/2
#categorize the new risk maps
brk  <- seq(0,1,0.2)
risk_hist_class <- raster::cut(risk_hist, breaks=brk) 
risk_ssp_class <- raster::cut(risk_ssp, breaks=brk) 

#overly FDSP locations on the new risk map 
PoC_points<-st_read('./input_data/spatial_data/wrl_prp_p_unhcr_PoC.geojson')
PoC_points$new_col<-c('PoC')
PoC_points$new_col<-factor(PoC_points$new_col)
#visualize the climate hazard class maps, as well as the new risk class maps with location of FDSP
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])

composite_hist_class_tbl<-overlayRast(composite_hist_class)
Visu_raster_discrete_with_points("./plots/with_FDSP/composite_hist_class_with_FDSP.tiff","",'Hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_hist_class_tbl,PoC_points)

composite_ssp_class_tbl<-overlayRast(composite_ssp_class)
Visu_raster_discrete_with_points("./plots/with_FDSP/composite_ssp_class_with_FDSP.tiff","",'Hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_ssp_class_tbl,PoC_points)

risk_hist_class_tbl<-overlayRast(risk_hist_class)
Visu_raster_discrete_with_points("./plots/with_FDSP/climate_gov_hist_with_FDSP.tiff","",'Risk class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),risk_hist_class_tbl,PoC_points)

risk_ssp_class_tbl<-overlayRast(risk_ssp_class)
Visu_raster_discrete_with_points("./plots/with_FDSP/climate_gov_ssp_with_FDSP.tiff","",'risk class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),risk_ssp_class_tbl,PoC_points)

#saving data of the new risk 
writeRaster(risk_hist_class, filename = ("./output_data/spatial_data/composite_gov_hist_class.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(risk_ssp_class, filename = ("./output_data/spatial_data/composite_gov_ssp_class.tif"), format = "GTiff",overwrite=TRUE)

#risk per country (and per class)
risk_data_hist<-raster::extract(risk_hist, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_hist<-risk_hist
values(empty_rast_hist)<-NA
the_risk_avrg<-risk_data_hist[dim(risk_data_hist)[2]]
risk_hist_avrg<-rasterize(the_risk_avrg, empty_rast_hist, field = the_risk_avrg@data[,1], update = TRUE)

risk_data_ssp<-raster::extract(risk_ssp, HoA,fun=mean, na.rm=TRUE,sp=TRUE) 
empty_rast_ssp<-risk_ssp
values(empty_rast_ssp)<-NA
the_risk_avrg2<-risk_data_ssp[dim(risk_data_ssp)[2]]
risk_ssp_avrg<-rasterize(the_risk_avrg2, empty_rast_ssp, field = the_risk_avrg2@data[,1], update = TRUE)

brk  <- seq(0,1,0.2)
risk_hist_class_avrg <- raster::cut(risk_hist_avrg, breaks=brk) 
risk_ssp_class_avrg <- raster::cut(risk_ssp_avrg, breaks=brk) 

#read the FDSP number per country
PoC<-read_csv("./input_data/FDSP_population.csv",skip=14)
our_PoC<-PoC[PoC$Year==2022,]
PoC_data<-data.frame(our_PoC)
PoC_data<-cbind(country_name=PoC_data[,4],country_code=PoC_data[,5],PoC_data[,c(6:10,12)])
PoC_data[PoC_data=='-']<-NA
PoC_data[,6]<-as.numeric(PoC_data[,6])
PoC_number<-rowSums(PoC_data[,3:8],na.rm=TRUE)
PoC_data<-data.frame(country_name=PoC_data[,1],country_code=PoC_data[,2],PoC_number)
PoC_data<-aggregate(PoC_number~(country_code+country_name),data=PoC_data,FUN=sum)
PoC_num<-NA
HoA<-cbind(HoA,PoC_num)
for (country_ind in (1:length(PoC_data$country_code))){
  the_ind<-which(HoA$iso3==(PoC_data$country_code[country_ind]))
  HoA$PoC_num[the_ind]<-PoC_data[country_ind,3]
  if (length(the_ind)==0){
    the_ind<-grep((PoC_data$country_name)[country_ind],HoA$gis_name)
    HoA$PoC_num[the_ind]<-PoC_data[country_ind,3]
  }
}

#build a combined shapefile with all regions
regional_files<-c('AMERICAs_15m_UNHCR.gpkg',"ASIA_15m_UNHCR.gpkg","East_and_horn_of_Africa_shapefile.gpkg",
                  "Europe_new.gpkg","MENA_polygon_2.gpkg","Southern_Africa_Ash.gpkg","West_and_central_new.gpkg")
all_region<-c('Americas','Asia','East and Horn of Africa','Europe','MENA','Southern Africa','West and Central Africa')
all_region_short<-c('Americas','Asia','EHAGL-Africa','Europe','MENA','S-Africa','WC-Africa')
all_region_short_bis<-c('Americas','Asia','EHAGL\nAfrica','Europe','MENA','S-Africa','WC\nAfrica')

for (count_region in 1:length(all_region)){
  polygon_rst <- st_read(paste0("./input_data/spatial_data/",regional_files[count_region]))
  the_region<-polygon_rst
  the_region$region_name<-rep(all_region[count_region],nrow(the_region))
  the_region$region_name_short<-rep(all_region_short[count_region],nrow(the_region))
  
  if (count_region==1)
  {combined_region<-the_region}
  if (count_region>1){
    combined_region<-rbind(combined_region,the_region)
  }
}

#add the region and hazard class and risk class for each country within the PoC data 
#########################
PoC_dict<-read_excel("./input_data/PoC_country_region_dictionary.xlsx")
PoC_dict<-as.data.frame(PoC_dict)

region<-c()
haz_class_hist<-c()
haz_class_ssp<-c()
haz_hist<-c()
haz_ssp<-c()
gov_value<-c()
risk_value_hist<-c()
risk_value_ssp<-c()
risk_class_hist<-c()
risk_class_ssp<-c()

for (i in 1:nrow(PoC_data)){
  
  the_ind<-which(HoA$iso3==(PoC_data$country_code[i]))
  if (length(the_ind)==0){
    the_ind<-grep((PoC_data$country_name)[i],HoA$gis_name)
  }
  if (length(the_ind)!=0){
    the_ind1<-the_ind[1]
    the_country_code<-(HoA$iso3)[the_ind1]
    the_reg<-(combined_region$region_name)[(combined_region$iso3)==the_country_code]
    if (length(the_reg)!=0){
      region<-c(region,the_reg)
    }else{
      #look in the other dictionary for region
      the_reg<-PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[i],2]
      region<-c(region,the_reg)
    }
    #hazard class per country
    all_haz_hist<-unlist(raster::extract(hist_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz1<- (as.integer(names(sort(table(all_haz_hist),decreasing=TRUE)[1])))
    if (length(the_haz1)==0){
      haz_class_hist<-c(haz_class_hist,NA)
    }else{
      haz_class_hist<-c(haz_class_hist,the_haz1)
    }
    
    all_haz_ssp<-unlist(raster::extract(ssp_class_avrg,HoA[the_ind,],exact=TRUE))
    the_haz2<-as.integer(names(sort(table(all_haz_ssp),decreasing=TRUE)[1]))
    if (length(the_haz2)==0){
      haz_class_ssp<-c(haz_class_ssp,NA)
    }else{
      haz_class_ssp<-c(haz_class_ssp,the_haz2)
    }
    #hazard average per country
    all_haz_hist_avr<-unlist(raster::extract(hist_avrg,HoA[the_ind,],exact=TRUE))
    the_haz1_avr<- (as.numeric(names(sort(table(all_haz_hist_avr),decreasing=TRUE)[1])))
    if (length(the_haz1_avr)==0){
      haz_hist<-c(haz_hist,NA)
    }else{
      haz_hist<-c(haz_hist,the_haz1_avr)
    }
    
    all_haz_ssp_avr<-unlist(raster::extract(ssp_avrg,HoA[the_ind,],exact=TRUE))
    the_haz1_avr<- (as.numeric(names(sort(table(all_haz_ssp_avr),decreasing=TRUE)[1])))
    if (length(the_haz1_avr)==0){
      haz_ssp<-c(haz_ssp,NA)
    }else{
      haz_ssp<-c(haz_ssp,the_haz1_avr)
    }
    
    #risk class per country
    all_risk_hist<-unlist(raster::extract(risk_hist_class_avrg,HoA[the_ind,],exact=TRUE))
    the_risk1<- (as.integer(names(sort(table(all_risk_hist),decreasing=TRUE)[1])))
    if (length(the_risk1)==0){
      risk_class_hist<-c(risk_class_hist,NA)
    }else{
      risk_class_hist<-c(risk_class_hist,the_risk1)
    }
    all_risk_ssp<-unlist(raster::extract(risk_ssp_class_avrg,HoA[the_ind,],exact=TRUE))
    the_risk2<- (as.integer(names(sort(table(all_risk_ssp),decreasing=TRUE)[1])))
    if (length(the_risk2)==0){
      risk_class_ssp<-c(risk_class_ssp,NA)
    }else{
      risk_class_ssp<-c(risk_class_ssp,the_risk2)
    }
    #risk average per country
    all_risk_hist_avr<-unlist(raster::extract(risk_hist_avrg,HoA[the_ind,],exact=TRUE))
    the_risk1_avr<- (as.numeric(names(sort(table(all_risk_hist_avr),decreasing=TRUE)[1])))
    if (length(the_risk1_avr)==0){
      risk_value_hist<-c(risk_value_hist,NA)
    }else{
      risk_value_hist<-c(risk_value_hist,the_risk1_avr)
    }
    all_risk_ssp_avr<-unlist(raster::extract(risk_ssp_avrg,HoA[the_ind,],exact=TRUE))
    the_risk2_avr<- (as.numeric(names(sort(table(all_risk_ssp_avr),decreasing=TRUE)[1])))
    if (length(the_risk2_avr)==0){
      risk_value_ssp<-c(risk_value_ssp,NA)
    }else{
      risk_value_ssp<-c(risk_value_ssp,the_risk2_avr)
    }

  }else{
    #look in the other dictionary for region
    the_reg<-PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[i],2]
    region<-c(region,the_reg)
    haz_class_hist<-c(haz_class_hist,NA)
    haz_class_ssp<-c(haz_class_ssp,NA)
    haz_hist<-c(haz_hist,NA)
    haz_ssp<-c(haz_ssp,NA)
    risk_class_hist<-c(risk_class_hist,NA)
    risk_class_ssp<-c(risk_class_ssp,NA)
    risk_value_hist<-c(risk_value_hist,NA)
    risk_value_ssp<-c(risk_value_ssp,NA)
  }

  the_ind_gov<-grep(PoC_data$country_code[i],gov_data$country_code)
  if (length(the_ind_gov)==0){
    the_ind_gov<-grep(PoC_data$country_name[i],gov_data$country_name_UNHCR)
  }
  if (length(the_ind_gov)!=0){
    gov_value<-c(gov_value,gov_data[the_ind_gov[1],5])
  }else{
    gov_value<-c(gov_value,NA)
  }
 
}
PoC_data<-cbind(PoC_data,region,haz_class_hist,haz_class_ssp,gov_value,haz_hist,haz_ssp,risk_class_hist,risk_class_ssp,risk_value_hist,risk_value_ssp)

saveRDS(PoC_data,'./output_data/FDSP_population_data_frame.Rds')

#Plot FDSP numbers per hazard class, overly each point with distribution over regions
region_col<-hcl.colors(length(all_region),palette='Dynamic')
names(region_col)<-all_region_short
the_titles<-c('Baseline','Future')
temporal<-c('hist','ssp')
for (t in (1:length(temporal))){
  class_count_ar<-c()
  for (cla in 1:5){

    all_PoC_class<-eval(parse(text=paste0('PoC_data[which(PoC_data$haz_class_',temporal[t],'==cla),]')))
    PoC_class_count<-sum(all_PoC_class$PoC_number,na.rm=TRUE)
    class_count_ar<-c(class_count_ar,PoC_class_count)
    
    #for each hazard class, plot a pie chart showing region
    reg_array<-c()
    for (the_reg in all_region){
      reg_array<-c(reg_array,sum(all_PoC_class$PoC_number[all_PoC_class$region==the_reg]))
    }
    the_df<-data.frame(all_region_short,reg_array)
    if (sum(reg_array)!=0){
      the_lab=paste0(round(reg_array / sum(reg_array) * 100, 1), "%")
    }else
    {
      the_lab=paste0(all_region_short, "0%")
    }
    the_df<-the_df[the_df$reg_array!=0 & the_lab!='0%',]
    the_lab<-the_lab[reg_array!=0 & the_lab!='0%']

    if (nrow(the_df)!=0){

      the_df2 <- the_df %>%
        mutate(csum = rev(cumsum(rev(reg_array))),
               pos = reg_array/2 + lead(csum, 1),
               pos = if_else(is.na(pos), reg_array/2, pos))

      piechart <- ggplot(the_df, aes(x="", y=reg_array, fill=as.factor(all_region_short))) +
        geom_col(color = "black") +
        coord_polar("y", start=0) +
        scale_fill_manual(values=region_col, name='region')+
        geom_label_repel(data = the_df2,
                         aes(y = pos, label = the_lab),
                         size = 13,nudge_x = 0.3,nudge_y = 0.9, show.legend = FALSE,force=5,segment.alpha = 0) +
        theme_void()+
        ggtitle(paste0('Hazard class:',cla))+
        theme(plot.title = element_text(hjust=0.5,size=15,face='bold'))+
        theme(legend.position = "none")
      tiff(paste0('./plots/with_FDSP/FDSP_number_',temporal[t],'_class_',toString(cla),'.tiff'), units="in", width=3, height=3.5, res=200,compression = 'lzw')
      print (piechart)
      dev.off()

    }

  }
  tiff(paste0('./plots/with_FDSP/FDSP_number_per_class_',temporal[t],'.tiff'),units="in", width=16, height=8, res=500,compression = 'lzw')
  par(mar=c(4,5.5,4,2) + 0.1)
  barplot(class_count_ar/1e6,width=1, names.arg = c('low','moderate','high','severe','extreme'), ylab="FDSP number (million)", xlab="Hazard class",  cex=2.2,xlim=c(0.3,5.7),
          ylim=c(0,max(class_count_ar/1e6)+25),cex.lab=3, cex.axis=2)
  dev.off()
  

}

#Plot number of FDSP per regions, overly each point with distribution over hazard classes
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
for (t in (1:length(temporal))){
  region_count_ar<-c()
  for (the_reg in 1:length(all_region)){
    reg<-all_region[the_reg]
    region_count_ar<-c(region_count_ar,sum(PoC_data$PoC_number[PoC_data$region==reg]))
    
    ##for each region, plot a pie chart showing hazard class
    the_reg_data<-PoC_data[PoC_data$region==reg,]

    class_array<-c()
    hazard_class<-1:5
    for (cla in hazard_class){
      region_ind<-eval(parse(text=paste0('the_reg_data[which(the_reg_data$haz_class_',temporal[t],'==cla),]')))
      region_ind<-sum(region_ind$PoC_number)
      class_array<-c(class_array,region_ind)
      
    }
    class_array<-round(class_array / sum(class_array) * 100, 1)
    class_array[which.max(class_array)]<- round(100-sum(class_array[-which.max(class_array)]),1)
    hazard_class<-hazard_class[class_array!=0]
    class_array<-class_array[class_array!=0]


    the_df<-data.frame(hazard_class,class_array)
    the_lab=paste0(as.character(the_df$class_array), "%")

    if (nrow(the_df)!=0){

      the_df2 <- the_df %>%
        mutate(csum = rev(cumsum(rev(class_array))),
               pos = class_array/2 + lead(csum, 1),
               pos = if_else(is.na(pos), class_array/2, pos))

      piechart <- ggplot(the_df, aes(x="", y=class_array, fill=as.factor(hazard_class))) +
        geom_col() +
        coord_polar("y", start=0) +
        scale_fill_manual(values=my_pals, name='hazard class')+
        geom_label_repel(data = the_df2,
                         aes(y = pos, label = the_lab),
                         size = 13, segment.alpha = 0) +
        theme_void()+
        ggtitle(all_region[the_reg])+
        theme(plot.title = element_text(hjust=0.5,size=15,face='bold'),panel.background = element_rect(fill = 'grey'))+
        theme(legend.position = "none")
      tiff(paste0('./plots/with_FDSP/FDSP_number_',temporal[t],'_',all_region[the_reg],'.tiff'), units="in", width=3, height=3.5, res=200,compression = 'lzw')

      print (piechart)
      dev.off()

    }

  }
  tiff(paste0('./plots/with_FDSP/FDSP_number_per_region_',temporal[t],'.tiff'),units="in", width=8, height=17, res=500,compression = 'lzw')
  par(mar=c(7,6,4,2) + 0.1)
  barplot(region_count_ar/1e6,width=0.97, names.arg = all_region_short_bis, xlab="", ylab="",  cex=3,ylim=c(0.3,7.5),
          xlim=c(0,max(region_count_ar/1e6)+5),cex.lab=3, cex.axis=2.2,horiz=TRUE)
  title(xlab = "FDSP number (million)", cex.lab = 3.5,
        line = 4.5)
  
  dev.off()
  
}
#############################################################################################
#Plot FDSP numbers per risk (including goverbance) class, overly each point with distribution over regions
region_col<-hcl.colors(length(all_region),palette='Dynamic')
names(region_col)<-all_region_short

for (t in (1:2)){

  class_count_ar<-c()
  for (cla in 1:5){
    all_PoC_class<-eval(parse(text=paste0('PoC_data[which(PoC_data$risk_class_',temporal[t],'==cla),]')))
    PoC_class_count<-sum(all_PoC_class$PoC_number,na.rm=TRUE)
    class_count_ar<-c(class_count_ar,PoC_class_count)
    
    #for each hazard class, plot a pie chart showing region
    reg_array<-c()
    for (the_reg in all_region){
      reg_array<-c(reg_array,sum(all_PoC_class$PoC_number[all_PoC_class$region==the_reg]))
    }
    the_df<-data.frame(all_region_short,reg_array)
    if (sum(reg_array)!=0){
      the_lab=paste0(round(reg_array / sum(reg_array) * 100, 1), "%")
    }else
    {
      the_lab="0%"
    }
    
    the_df<-the_df[the_df$reg_array!=0 & the_lab!='0%',]
    the_lab<-the_lab[reg_array!=0 & the_lab!='0%']
    
    if (nrow(the_df)!=0){
      
      the_df2 <- the_df %>%
        mutate(csum = rev(cumsum(rev(reg_array))),
               pos = reg_array/2 + lead(csum, 1),
               pos = if_else(is.na(pos), reg_array/2, pos))
      
      piechart <- ggplot(the_df, aes(x="", y=reg_array, fill=as.factor(all_region_short))) +
        geom_col(color = "black") +
        coord_polar("y", start=0) +
        scale_fill_manual(values=region_col, name='region')+
        geom_label_repel(data = the_df2,
                         aes(y = pos, label = the_lab),
                         size = 13,nudge_x = 0,nudge_y = 0, show.legend = TRUE,force=5,segment.alpha = 0,max.overlaps=20) +
        theme_void()+
        ggtitle(paste0('risk class:',cla))+
        theme(plot.title = element_text(hjust=0.5,size=15,face='bold'))+
        theme(legend.position = "none")
      tiff(paste0('./plots/with_FDSP/FDSP_number_gov_',temporal[t],'_class_',toString(cla),'.tiff'), units="in", width=3, height=3.5, res=200,compression = 'lzw')
      print (piechart)
      dev.off()
      
    }
    
  }
  

  tiff(paste0('./plots/with_FDSP/FDSP_number_gov_per_class_',temporal[t],'.tiff'),units="in", width=16, height=8, res=500,compression = 'lzw')
  par(mar=c(4,5.5,4,2) + 0.1)
  barplot(class_count_ar/1e6,width=1, names.arg = c('low','moderate','high','severe','extreme'), ylab="FDSP number (million)", xlab="Hazard class",  cex=2.2,xlim=c(0.3,5.7),
          ylim=c(0,max(class_count_ar/1e6)+25),cex.lab=3, cex.axis=2)
  dev.off()
}


#Plot number of FDSP per regions, overly each point with distribution over risk classes (with governance)
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
for (t in (1:length(temporal))){
  region_count_ar<-c()
  for (the_reg in 1:length(all_region)){
    reg<-all_region[the_reg]
    region_count_ar<-eval(parse(text=paste0('c(region_count_ar,sum(PoC_data$PoC_number[which(PoC_data$region==reg & !is.na(PoC_data$risk_class_',temporal[t],') )]))')))
    
    ##for each region, plot a pie chart showing hazard class
    the_reg_data<-eval(parse(text=paste0('PoC_data[which(PoC_data$region==reg & !is.na(PoC_data$risk_class_',temporal[t],') ),]')))
    class_array<-c()
    the_class<-1:5
    for (cla in the_class){
      region_ind<-eval(parse(text=paste0('the_reg_data[which(the_reg_data$risk_class_',temporal[t],'==cla),]')))
      region_ind<-sum(region_ind$PoC_number)
      class_array<-c(class_array,region_ind)
    }
    class_array<-round(class_array / sum(class_array) * 100, 1)
    class_array[which.max(class_array)]<- round(100-sum(class_array[-which.max(class_array)]),1)
    the_class<-the_class[class_array!=0]
    class_array<-class_array[class_array!=0]
  
    the_df<-data.frame(the_class,class_array)
    the_lab=paste0(as.character(the_df$class_array), "%")
    
    if (nrow(the_df)!=0){
      
      the_df2 <- the_df %>%
        mutate(csum = rev(cumsum(rev(class_array))),
               pos = class_array/2 + lead(csum, 1),
               pos = if_else(is.na(pos), class_array/2, pos))
      
      piechart <- ggplot(the_df, aes(x="", y=class_array, fill=as.factor(the_class))) +
        geom_col() +
        coord_polar("y", start=0) +
        scale_fill_manual(values=my_pals, name='risk class')+
        geom_label_repel(data = the_df2,
                         aes(y = pos, label = the_lab),
                         size = 13, segment.alpha = 0) +
        theme_void()+
        ggtitle(all_region[the_reg])+
        theme(plot.title = element_text(hjust=0.5,size=15,face='bold'))+
        theme(legend.position = "none")
      tiff(paste0('./plots/with_FDSP/FDSP_number_gov_',temporal[t],'_',all_region[the_reg],'.tiff'), units="in", width=3, height=3.5, res=200,compression = 'lzw')
      
      print (piechart)
      dev.off()
      
    }
  }
  tiff(paste0('./plots/FDSP_number_gov_per_region ',temporal[t],'.tiff'),units="in", width=8, height=17, res=500,compression = 'lzw')
  par(mar=c(7,6,4,2) + 0.1)
  barplot(region_count_ar/1e6,width=0.93, names.arg = all_region_short_bis, xlab="", ylab="",  cex=3,ylim=c(0.3,7.5),
          xlim=c(0,max(region_count_ar/1e6)+4),cex.lab=3, cex.axis=2.2, horiz=TRUE)
  title(xlab = "FDSP number (million)", cex.lab = 3.5,
        line = 4.5)
  dev.off()
   
}
