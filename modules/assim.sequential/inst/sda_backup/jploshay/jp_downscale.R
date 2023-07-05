library(terra)
library(sp)
library(randomForest)
library(tidyverse)
library(sf)
library(ggmap)
library(hydroGOF)

### Setup input data

# SpatRater stack of covariate maps (average temperature, total precip. solar radiation, and vapor pressure)
clim_files <- rast(c("/projectnb/dietzelab/jploshay/pecan_copy/jploshay/avg_temp_prep/tavg_10m.tif",
                     "/projectnb/dietzelab/jploshay/pecan_copy/jploshay/total_prec/precip_10m.tif",
                     "/projectnb/dietzelab/jploshay/pecan_copy/jploshay/10m_srad/10m_srad_map.tif",
                     "/projectnb/dietzelab/jploshay/pecan_copy/jploshay/10m_vapr/10m_vapr_map.tif"))
names(clim_files) <- c("tavg", "prec", "srad", "vapr")


# Pull site identification numbers 
analysis <- readRDS("/projectnb/dietzelab/jploshay/pecan_copy/jploshay/downscaling/data/ANALYSIS.rds")
sites <- unique(attributes(analysis[[names(analysis)[1]]])$Site)


# Pull SDA ID, site name, longitude, and latitude 
site_info <- read_csv("/projectnb/dietzelab/jploshay/pecan_copy/jploshay/downscaling/data/sda_500_flux_net_site_table.csv")
site_info <- site_info %>%
  select(sda_id, lon, lat, site_name)
site_info$sda_id <- as.character(site_info$sda_id)


# Isolate above ground wood data from 2012 into a dataframe
index = which(names(analysis) == "2012/07/15") 
data = analysis[[index]]
data = data[which(names(data) == "AbvGrndWood")]


# Create data frame of the median above ground wood data binded to the SDA site ID
biomass = apply(data, 2, FUN = median, na.rm = T)
biomass = as.data.frame(cbind(sites, biomass))
row.names(biomass) = NULL
names(biomass) = c("site_id", "biomass")



biomass_data <- inner_join(biomass, site_info, by = c("site_id" = "sda_id")) 
coordinates(biomass_data) = ~lon+lat
sp::proj4string(biomass_data) <- sp::CRS("+init=epsg:4326")
biomass_data <- sp::spTransform(biomass_data, crs(clim_files))#reproject


biomass_for_plot <- sf::st_as_sf(biomass_data)
biomass_for_model <- as.data.frame(biomass_for_plot)

biomass_map <- get_stamenmap(
  bbox = unname(st_bbox(biomass_for_plot)),
  zoom = 5, maptype = 'toner-lite', source = 'stamen'
) %>% ggmap()

biomass_map + 
  geom_sf(
    data = biomass_for_plot, 
    aes(size = 1), 
    color = 'red', alpha = 0.5,
    show.legend = 'point', inherit.aes = F
  )

# prepare training and testing datasets
# lon_lat <- cbind(analysis_data$lon, analysis_data$lat)
# lon_lat <- vect(lon_lat, crs = sf::st_crs(4326)[[2]])
# 
# t = as.data.frame(terra::extract(clim_files, lon_lat))


# bionas_vect <- vect(biomass_for_plot, geom="wkt", crs="+init=epsg:4326", keepgeom=F)
# input_data = as.data.frame(extract(clim_files, biomass_for_plot)) #???
# input_data = cbind(input_data, biomass$biomass)

lon_lat <- cbind(site_info$lon, site_info$lat)
lon_lat <- vect(lon_lat, crs = sf::st_crs(4326)[[2]])

input_data = as.data.frame(terra::extract(clim_files, lon_lat))
input_data = cbind(input_data, biomass_data$biomass)
names(input_data) = c("site", "tavg", "prec", "srad", "vapr", "biomass")
input_data$biomass = as.numeric(input_data$biomass)


#### outlier analysis ####
# sd = sd(input_data$biomass)
# z = abs(input_data$biomass-mean(input_data$biomass))/sd
# index = which(z > 2.5)
# input_data = input_data[-index,]


samples = sample(1:nrow(input_data), .8*nrow(input_data))
training = input_data[samples,]
testing = input_data[-samples,]


# make the model!!!
rf = randomForest(biomass~tavg+prec+srad+vapr, 
                  data = training, ntree = 1000, na.action = na.omit, keep.forest = T)

pr = predict(rf, newdata = testing, na.action = na.pass)
max = max(pr, testing$biomass)
# jpeg('/data2/bmorrison/sda/500_site_run/rf_model_diagnostics_2021.jpeg', height = 8, width = 8, units = "in", res = 300)
par(mfrow = c(2,2))
plot(rf, main = "Tree Error")
varImpPlot(rf, main = "Variable Importance")
plot(testing$biomass, pr, xlim = c(0, max), ylim = c(0, max), xlab = "SDA AGB", ylab = "Predicted", main = "Obs vs. Predicted")
abline(0, 1, col = 'red')

c = round(cor(pr, testing$biomass, use ="complete.obs"), digits = 1)
rmse = function(obs, pred)
{
  error = sqrt(sum((pred-obs)^2)/length(pred))
  return(error)
}
error = round(rmse(testing$biomass, pr), digits = 1)
pb = pbias(testing$biomass, pr, na.rm = T)
text(x = 0, y = 200, labels = paste0("Cor=",c), pos = 4)
text(x = 0, y = 190, labels = paste0("RMSE=", error), pos = 4)
text(x = 0, y = 180, labels = paste0("%Bias=",pb), pos = 4)

test = predict(object = clim_files, model = rf, na.rm = T)
plot(test, main = "CONUS AGB Estimate")
dev.off()

##### Make STDEV estimate from model
na_index = which(is.na(clim_files[]))
d = as.data.frame(clim_files)
d$cell = 1:nrow(d)
d = d[-na_index,]
test2 = predict(rf, newdata = d, na.action = na.omit, predict.all = T)
test2_data = test2$individual
test2_sd = apply(test2_data, 1, FUN = sd, na.rm = T)

test2 = mask*NA
test2[d$cell] = test2_sd

lt_agb = raster("/data2/bmorrison/sda/500_site_run/data/landtrendr_agb_2012_800m.tif")
lt_agb_se = raster("/data2/bmorrison/sda/500_site_run/data/landtrendr_agb_stdev_2012_800m.tif")
agb_diff = abs(lt_agb-test)
se_diff = abs(lt_agb_se-test2)
# range_lt = range(r[], na.rm = T)
# breaks_lt = seq(range_lt[1], range_lt[2], by = 31)
breaks_agb = c(0, 50, 100, 150 200, 250, 300, 350, 400, 450, 500, 550, 600)
breaks_stdev = c(0, 25, 50, 75, 100, 125, 150)
breaks_diff = c()
jpeg('/data2/bmorrison/sda/500_site_run/rf_model_comparison_2021.jpeg', height = 10, width = 14, units = "in", res = 300)
par(mfrow = c(2,3))
plot(lt_agb, col = rev(terrain.colors(length(breaks_agb)-1)), breaks = breaks_agb, main = "LandTrendr AGB 800m")
plot(lt_agb_se, col = rev(terrain.colors(length(breaks_stdev)-1)), breaks = breaks_stdev, main = "LandTrendr STDEV 800m")
plot(agb_diff,  main = "Difference LT vs. RF AGB Estimates")
plot(test, col = rev(terrain.colors(length(breaks_agb)-1)), breaks = breaks_agb, main = "SDA RF AGB")
plot(test2, col = rev(terrain.colors(length(breaks_stdev)-1)), breaks = breaks_stdev, main = "SDA RF STDEV")
plot(se_diff,  main = "Difference LT vs. RF STDEV Estimates")
dev.off()