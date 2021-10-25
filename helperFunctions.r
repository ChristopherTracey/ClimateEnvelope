# Here's a function to plot the current and future HSM projections
project.sdm <- function(prediction, plotName){
  rng = range(c(0, 1)) #a range to have the same min and max for both plots
  map.p <- rasterToPoints(prediction)
  df <- data.frame(map.p) # Make the points a dataframe for ggplot
  colnames(df) <- c('Longitude', 'Latitude', 'Probability') # Make appropriate column headings
  ggplot() +
    geom_raster(data=df,  mapping=aes(y=Latitude, x=Longitude, fill=Probability)) +
    scale_fill_gradient2(low="blue", high="red", guide="colorbar", limits=c(floor(rng[1]), ceiling(rng[2]))) +
    geom_sf(data=studyArea, color='black', fill=NA) +
    geom_point(data=sp_coords, aes(x=lon, y=lat), color='gray', size=2, shape=4) +
    ggtitle(plotName) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()  
    )
  
}
