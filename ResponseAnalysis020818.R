require(openxlsx)
require(raster)
require(downloader)
require(rgdal)
require(sp)
require(raster)
require(quantreg)
require(leaflet)
require(htmlwidgets)
require(dplyr)
require(randomForestSRC)
require(lme4)


##add rainfall (optional)
# rain.data<-read.csv("data/processed/CHIRPSdata.csv")
# surveydata<-cbind(surveydata, rain.data)
# ##assign monthly rainfall for each year
# monthly.rain<-matrix(,nrow(surveydata),12)
# colnames(monthly.rain)<-paste("rain",1:12,sep="_")
# sel<-intersect(grep("rain",colnames(surveydata)),grep("2016",colnames(surveydata))) ##columns with 2016 rain
# monthly.rain[which(surveydata$Year=="2016"),]<-as.matrix(surveydata[which(surveydata$Year=="2016"), sel])
# sel<-intersect(grep("rain",colnames(surveydata)),grep("2017",colnames(surveydata))) ##columns with 2016 rain
# monthly.rain[which(surveydata$Year=="2017"),]<-as.matrix(surveydata[which(surveydata$Year=="2017"), sel])
# surveydata<-surveydata[,grep("rain_",colnames(surveydata),invert=T)]
# surveydata<-cbind(surveydata, monthly.rain)
# Create a data folder in your current working directory
dir.create("OAF_data", showWarnings=F)
setwd("./OAF_data")


##get R object with yield data
download("https://www.dropbox.com/s/yuaceydd6fg1auf/clean_yga_val_data.Rdata?raw=1", "clean_yga_val_data.Rdata", mode = "wb")

load("clean_yga_val_data.Rdata")


# download GADM-L3 shapefile (@ http://www.gadm.org)
download("https://www.dropbox.com/s/otspr9b9jtuyneh/KEN_adm3.zip?raw=1", "KEN_adm3.zip", mode = "wb")
unzip("KEN_adm3.zip", overwrite = T)
shape <- shapefile("KEN_adm3.zip")

# download raster stack
download("https://www.dropbox.com/s/mz1t0zyq8uoqrhq/KE_250m_2017.zip?raw=1", "KE_250m_2017.zip", mode = "wb")
unzip("KE_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)


#set data
surveydata <- surveydata %>% dplyr::select(-one_of(c('BoxALength', 'BoxAWidth', 'BoxBLength', 'BoxBWidth', 'DryWeightBagA', 'DryWeightBagB', 'acc', 'CMID')))
columns <-c("alt", "CANKGS", "DAPKGS", "PlotSize", "KeepCropSeparate", "AddedOrRemovedCropFromBag", 'DAPMaize','lon', 'lat')
surveydata[, columns] <- lapply(columns, function(x) as.numeric(as.character(surveydata[[x]])))
# add N and P totals to data 

surveydata <- surveydata %>% mutate(
  tot_N=(0.18*DAPKGS)+(0.26*CANKGS), 
  tot_P=(0.46*DAPKGS)*0.44, 
  dap.kg.per.acre = tot_N/PlotSize, 
  total_N.per.ha = dap.kg.per.acre*2.47105, 
  can.kg.per.acre = tot_P/PlotSize, 
  total_P.per.ha = can.kg.per.acre*2.47105) 

yield.data<-surveydata

sel.cols<-unique(c(which(colnames(yield.data)%in%c("yield","Year","PlotSize","total_N.per.ha","total_P.per.ha","lon","lat","alt")),grep("rain_",colnames(yield.data))))
yield.data<-yield.data[,sel.cols]
#remove missing data
yield.data<-na.omit(yield.data)
yield.data<-unique(yield.data)
##fix lat long
yield.data<-plyr::rename(yield.data, c("lat" = "lon","lon"="lat"))
yield <- yield.data


# set ROI grid extent
ext <- data.frame(lat = c(-1.2,-1.2,1.2,1.2), lon = c(33.9,35.5,33.9,35.5)) ## set ROI extent in degrees
names(ext) <- c("lat","lon")
coordinates(ext) <- ~ lon + lat
proj4string(ext) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ext <- spTransform(ext, CRS("+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
bb <- extent(ext)
grids <- crop(grids, bb)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(yield) <- ~lon+lat
projection(yield) <- projection(shape)
gadm <- yield %over% shape
yield <- as.data.frame(yield)
yield <- cbind(gadm[ ,c(5,7,9)], yield)
#colnames(yield) <- c("district","division","location","id","lat","lon","year","trt","can","dap","fsize","yield")

# project survey coords to grid CRS
yield.proj <- as.data.frame(project(cbind(yield$lon, yield$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(yield.proj) <- c("x","y")
yield <- cbind(yield, yield.proj)
coordinates(yield) <- ~x+y
projection(yield) <- projection(yield)

# extract gridded variables at survey locations
yieldgrid <- extract(grids, yield)
gsdat <- as.data.frame(cbind(yield, yieldgrid)) 
#remove x, y
gsdat$lat<-NULL
gsdat$lon<-NULL
gsdat<-na.omit(gsdat)

#####make frame with yield,year and predictors
total.dat<-gsdat[,-c(1:3)]


##subset based on criteria
total.dat<-total.dat[is.finite(total.dat$total_N.per.ha),]
total.dat<-total.dat[total.dat$yield<10,]
###decide if rain variables enter
#total.dat<-total.dat[,grep("rain_",colnames(total.dat),invert=T)]
###decide if non-spatial variables enter
#total.dat<-total.dat[,!colnames(total.dat)%in%c("PlotSize","total_N.per.ha","total_P.per.ha","alt")]



###make first attempt at prediction
test=0
pix.size=100
while(test==0){
#set pixel size to assure at least 3 N rates per pixel
total.dat$x.p<-round(total.dat$x/pix.size)*pix.size
total.dat$y.p<-round(total.dat$y/pix.size)*pix.size
total.dat$yxy.p<-paste(total.dat$Year,total.dat$x.p,total.dat$y.p)
tap<-tapply(total.dat$total_N.per.ha,total.dat$yxy.p,function(x) length(unique(x)))
tab<-table(tap)
tab<-tab/sum(tab)
test<-tab[1]<0.2	
pix.size=pix.size+10
}

##remove site/years with less than 2 N levels
sel.sites<-tap[tap>1] #select pixels with more than one N value
total.dat<-total.dat[total.dat$yxy.p%in%names(sel.sites),]





###now estimate location specific response


dat2016<-total.dat[total.dat$Year=="2016",]
dat2017<-total.dat[total.dat$Year=="2017",]


##2016

dat2016$ID<-factor(dat2016$yxy.p)
dat2016$tot.N.scaled<-dat2016$total_N.per.ha/100 #scale N for convergence

###fit model and extract random slopes
#lmm16<-lmer(yield~ tot.N.scaled +total_P.per.ha+(tot.N.scaled ||ID),data= dat2016)
lmm16<-lmer(yield~ tot.N.scaled +(tot.N.scaled ||ID),data= dat2016)
response.frame<-ranef(lmm16)$ID
response.frame<-t(t(response.frame)+fixef(lmm16)[1:2]) ##add fixed effects
response.frame<-as.data.frame(response.frame)
colnames(response.frame)<-c("base.yield.est","response.est")
response.frame$ID<-factor(rownames(response.frame))
dat2016.means<-aggregate(dat2016,list(dat2016$ID),mean,na.rm=T)
dat2016.means$ID<-dat2016.means$Group.1
##merge with estimated responses
dat2016.means<-merge(dat2016.means, response.frame)
rownames(dat2016.means)<-dat2016.means$ID

##2017

dat2017$ID<-factor(dat2017$yxy.p)
dat2017$tot.N.scaled<-dat2017$total_N.per.ha/100 #scale N for convergence

###fit model and extract random slopes
#lmm17<-lmer(yield~ tot.N.scaled +total_P.per.ha+(tot.N.scaled ||ID),data= dat2017)
lmm17<-lmer(yield~ tot.N.scaled +(tot.N.scaled ||ID),data= dat2017)
response.frame<-ranef(lmm17)$ID
response.frame<-t(t(response.frame)+fixef(lmm17)[1:2]) ##add fixed effects
response.frame<-as.data.frame(response.frame)
colnames(response.frame)<-c("base.yield.est","response.est")
response.frame$ID<-factor(rownames(response.frame))
dat2017.means<-aggregate(dat2017,list(dat2017$ID),mean,na.rm=T)
dat2017.means$ID<-dat2017.means$Group.1
##merge with estimated responses
dat2017.means<-merge(dat2017.means, response.frame)
rownames(dat2017.means)<-dat2017.means$ID


#remove variables that should not be used in prediction
rm.var<-c("ID","Group.1","Year","alt","PlotSize","yield","total_N.per.ha","total_P.per.ha","x.p","y.p","yxy.p","tot.N.scaled","base.yield.est")

##extract district info
coordinates(dat2016.means) <- ~x+y
coordinates(dat2017.means) <- ~x+y
shape<-spTransform(shape,CRS("+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
projection(dat2016.means) <- projection(shape)
projection(dat2017.means) <- projection(shape)
gadm2016 <- dat2016.means %over% shape
gadm2017 <- dat2017.means %over% shape

##reconvert to normal data.frame
dat2016.means<-cbind(dat2016.means@data, coordinates(dat2016.means))
dat2017.means<-cbind(dat2017.means@data, coordinates(dat2017.means))


dat2016.means<-dat2016.means[,!colnames(dat2016.means)%in% rm.var]
dat2017.means<-dat2017.means[,!colnames(dat2017.means)%in% rm.var]


##fit RF model
rf.reg16 <- rfsrc(response.est ~., data = dat2016.means,importance = T)
cor(rf.reg16$yvar, rf.reg16$predicted.oob)^2


##fit RF model
rf.reg17 <- rfsrc(response.est ~., data = dat2017.means,importance = T)

##oob r-squared
cor(rf.reg17$yvar, rf.reg17$predicted.oob)^2

##check variable importance
rf.reg16$importance[order(rf.reg16$importance,decreasing=T)]
rf.reg17$importance[order(rf.reg17$importance,decreasing=T)]

##check correlation variable importance across years
#r-seuqred
cor(rf.reg17$importance, rf.reg16$importance)^2
#plots
plot(rf.reg17$importance, rf.reg16$importance)
plot(rank(rf.reg17$importance),rank(rf.reg16$importance))




#####predict across years

pred1<-predict(rf.reg16,newdata= dat2017.means)
pred2<-predict(rf.reg17,newdata= dat2016.means)

##r-squared plot level predictions/observations across years
cor(dat2017.means$response.est,pred1$predicted)^2
cor(dat2016.means$response.est,pred2$predicted)^2

##plots
plot(dat2017.means$response.est,pred1$predicted)
plot(dat2016.means$response.est,pred2$predicted)



##check accuracy when averaged over districts

tap.obs2017<-tapply(dat2017.means$response.est, gadm2017$NAME_1,mean)
tap.pred2017<-tapply(pred1$predicted, gadm2017$NAME_1,mean)

tap.obs2016<-tapply(dat2016.means$response.est, gadm2016$NAME_1,mean)
tap.pred2016<-tapply(pred2$predicted, gadm2016$NAME_1,mean)


##r-squared district level predictions/observations
cor(tap.obs2017, tap.pred2017)^2
cor(tap.obs2016, tap.pred2016)^2

#plot predictions
par(mfrow=c(2,2))
lm1<-lm(dat2017.means$response.est ~pred1$predicted)
plot(formula(lm1),xlab="response 2017, predicted from 2016",ylab="response 2017, observed",main="plot level")
abline(coef(lm1),col="blue")
lm2<-lm(dat2016.means$response.est ~pred2$predicted)
plot(formula(lm2),xlab="response 2016, predicted from 2017",ylab="response 2016, observed",main="plot level")
abline(coef(lm2),col="blue")
lm3<-lm(tap.obs2017~tap.pred2017)
plot(formula(lm3),pch=19,xlab="response 2017, predicted from 2016",ylab="response 2017, observed",main="district level")
abline(coef(lm3),col="blue")
lm4<-lm(tap.obs2016~tap.pred2016)
plot(formula(lm4),pch=19,xlab="response, predicted from 2017",ylab="response 2016, observed",main="district level")
abline(coef(lm4),col="blue")


###make gridded prediction surface for 2016
##make prediction on grid for both models and plot (for spatial variables only)
grid.dat<-na.omit(as.data.frame(grids,xy=T))
grid.dat.2016<-grid.dat.2017<-grid.dat
grid.dat.2016$Year<-2016
grid.dat.2017$Year<-2017


###predict
grid.predict.2017<-predict(rf.reg16,newdata= grid.dat.2017)
grid.predict.2016<-predict(rf.reg17,newdata= grid.dat.2016)

##make rasters for plotting
pred.rast.2017<-cbind(grid.predict.2017$xvar[,c("x","y")],grid.predict.2017$predicted)
colnames(pred.rast.2017)<-c("x","y","z")
pred.rast.2017<-rasterFromXYZ(pred.rast.2017)
pred.rast.2016<-cbind(grid.predict.2016$xvar[,c("x","y")],grid.predict.2016$predicted)
colnames(pred.rast.2016)<-c("x","y","z")
pred.rast.2016<-rasterFromXYZ(pred.rast.2016)


par(mfrow=c(1,2))
image(pred.rast.2017,main="2017 response, predicted from 2016")
points(dat2017$x,dat2017$y,col=gray(0.7, alpha = 0.1))
image(pred.rast.2016,main="2016 response, predicted from 2017")
points(dat2016$x,dat2016$y,col=gray(0.7, alpha = 0.1))

