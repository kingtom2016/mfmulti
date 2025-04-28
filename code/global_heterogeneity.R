library(raster)
library(terra)
library(sp)

source("D:/Myfile/Research/Script/load_package.txt")
# clim <- getData('worldclim', var='bio', res=10) # BIO1 for MAT
clim <- geodata::worldclim_global(var='bio', res=10, version="2.1",path="./")
data_gis<- rast("H:/analysis/test/新建文件夹/Entropy_01_05_1km_uint16.tif")# as.data.frame(data_gis) %>% head()
col=c("agricultural land"="#c49c94","broadleaf"="#00a14e","coniferous"="#FFA500","grassland"="#42b540","shrubland"="#767676","bare land"="#f8c9de","forest"="#ecc37a")
df<-df_bac
  points <- cbind(df$lon, df$lat)
  df$raster_value <- extract(data_gis,points)[,1]  #data_gis_euro
  df$MAT <- extract(clim,points)[,1]
  df$MAP <- extract(clim,points)[,12]
df<-df %>% filter(!is.na(value))
df<-df %>% mutate(raster_value=scale(raster_value))
res<-lmer(value ~ raster_value + MAT + MAP+ (1|subtype)  + (1|study_id) ,  data = df)
Anova(res)
sjPlot::tab_model(res,show.intercept = F,digits = 6,file = "lmer_bac.html")
df$value_preidct <- predict(res, df)
p1<-df %>% mutate(subtype=factor(subtype,levels=names( col))) %>% 
  ggplot(aes(x=raster_value,y=value)) + geom_point(aes(color=subtype),alpha=0.5,show.legend = TRUE) + 
  scale_color_manual(values=col,drop = F,guide=guide_legend(override.aes = list(size=2.5)))+
  geom_smooth(aes(color=subtype,y=value_preidct),method = "lm",show.legend = F)+
  scale_x_continuous(expand = expansion(c(0.15,0.05)))+
  labs(x="Normalized habitat heterogeneity",y="Shannon index",color="")+labs(title="Bacteria")+theme_minimal();p1  #+theme(legend.position = c(0.15,0.15))+guides(color=guide_legend(override.aes = list(alpha=0.5))) ;p1





df_fun<-p_fungi_shannon_entropy$data%>% filter(!is.na(value));df<-df_fun
points <- cbind(df$lon, df$lat)
df$MAT <- extract(clim,points)[,1]
df$MAP <- extract(clim,points)[,12]
df<-df %>% filter(!is.na(MAP))
df<-df %>% mutate(raster_value=scale(raster_value))
res<-lmer(value ~ raster_value + MAT + MAP +  (1|subtype) ,  data = df);length(predict(res, re.form = NULL)  )
Anova(res)
df$value_preidct <- predict(res, df)
sjPlot::tab_model(res,show.intercept = F,digits = 6,file = "lmer_fun.html")

p2<-df %>%  
  ggplot(aes(x=raster_value,y=value)) + geom_point(aes(color=subtype),alpha=0.5) + 
  scale_color_manual(values=col,drop = FALSE)+
  geom_smooth(aes(color=subtype,y=value_preidct),method = "lm")+
  labs(x="Normalized habitat heterogeneity",y="Shannon index",color="")+labs(title="Fungi")+
  theme_minimal()#+theme(legend.position = "none")
p2



df_meta<-p_meta_shannon_entropy$data ;df<-df_meta
points <- cbind(as.numeric(df$lon), as.numeric(df$lat))
df$MAT <- extract(clim,points)[,1]
df$MAP <- extract(clim,points)[,12]
df<-df %>% filter(!is.na(MAP)) 
df<-df %>% mutate(raster_value=scale(raster_value))
res<-lmer(value ~ raster_value + MAT + MAP + (1|Bioproject) + (1|subtype) , 
          data = df  );length(predict(res, re.form = NULL)  )
Anova(res)
df$value_preidct <- predict(res, df)
sjPlot::tab_model(res,show.intercept = F,digits = 6,file = "lmer_function.html")

p3<-df %>% 
  ggplot(aes(x=raster_value,y=value)) + geom_point(aes(color=subtype),alpha=0.5) + 
  scale_color_manual(values=col,drop = FALSE)+
  geom_smooth(aes(color=subtype,y=value_preidct),method = "lm")+
  labs(x="Normalized habitat heterogeneity",y="Metagenomic diversity index",color="")+labs(title="Metagenome")+
  theme_minimal()+coord_cartesian(ylim=c(0.5,1))#+theme(legend.position = "none");
p3


# p1+p2+p3+plot_layout(guides = "collect")
ggarrange(p1,p2,p3,nrow=1,common.legend = T,legend = "right" )
		  
		  
library(terra)
data_gis<- rast("H:/analysis/test/新建文件夹/Entropy_01_05_1km_uint16.tif")
colors <- colorRampPalette( colors = c("#85c1c8", "#af4980", "#c0182a", "#d33300", "#e99900",  "#ffff00"),bias = 0.25)
colors <- colorRampPalette( colors = c("white","#f9d1c1", "#fccdbb", "#f56954", "#de2826", "#8b191b"),bias = 0.25)
colors <- colorRampPalette( colors = c("white","#ebf9fe","#d9f4fe", "#a9e6fc", "#6dd5fa", "#57b9e5","#2980b9"),bias = 0.25)

pdf("p_global_heterogenity1.pdf",width=10,height=6)
par(mfrow=c(1,1))
plot(data_gis,col=colors(100),breaks = c(0, 5000, 10000, 15000,20000,25000,35000), legend="bottomleft", axes=FALSE, box=FALSE,title="Habitat heterogenity")
# title("Bacteria", adj=0.075, font.main = 1, cex.main = 1.5)
tmp<-df_fun %>% dplyr::select(lon,lat) %>% as.matrix()
pts <- vect(as.data.frame(tmp), geom=c("lon", "lat"), crs="+proj=longlat +datum=WGS84")
points(pts, col="#ffc312", cex=1, pch=20)
# dev.off()
# pdf("p_global_heterogenity2.pdf",width=12,height=3.5)
# plot(data_gis,col=colors(100),breaks = c(0, 5000, 10000, 15000,20000,25000,35000), legend=NA )
# title("Fungi", adj=0.075, font.main = 1, cex.main = 1.5)
tmp<-df_bac %>% dplyr::select(lon,lat) %>% as.matrix()
pts <- vect(as.data.frame(tmp), geom=c("lon", "lat"), crs="+proj=longlat +datum=WGS84")
points(pts, col="#e3adb9", cex=1, pch=20)
# pdf("p_global_heterogenity2.pdf",width=12,height=3.5)
# plot(data_gis,col=colors(100),breaks = c(0, 5000, 10000, 15000,20000,25000,35000), legend=NA )
# title("Function", adj=0.075, font.main = 1, cex.main = 1.5)
tmp<-df_meta %>% dplyr::select(lon,lat) %>% as.matrix()
pts <- vect(as.data.frame(tmp), geom=c("lon", "lat"), crs="+proj=longlat +datum=WGS84")
points(pts, col="#954024", cex=1, pch=20)
legend(x=-180,y=15, legend=c("Bacterial amplicon","Fungal amplicon", "Metagenome"), col=c("#ffc312", "#e3adb9", "#954024"), pch=20, cex=0.8,bty = "n", title="Datasets")
dev.off()