---
title: "Spatial Interpolation with R"
output: 
    html_document : default
---
<style>
pre{
    color: black;
}

pre:not([class]){
    color: white;
    background-color: #444;
}

</style>



```{r setup, include=FALSE}
library(knitr)
knitr::opts_chunk$set(echo = TRUE)

```
<p>Author: Liam Osler, Date: March 2021</p>
<div id="target"></div>

<p><b>Project Questions</b>,</p>
<a href = "#q1">Identify which interpolation method best approximated the continuous elevation surface and justify your answer.</a>
<a href = "#q2">Identify which interpolation method least approximated the continuous elevation surface and justify your answer.</a>


<h2>Initial setup:</h2>
<h3>
<a href="interpolation.rmd" download="">Download the .rmd notebook for this project
</a>
</h3>
<p>Install the required packages and their dependencies from the <a href="https://mirror.its.dal.ca/cran/">CRAN - The Comprehensive R Archive Network</a> if needed (manual documentation is also available there):</p>
```{r}
#Uncomment the following line to include the installation of project packages if desired:

#install.packages(c("dplyr", "leaflet","ggplot2", "ggrepel", "ggspatial", "libwgeom", "lwgeom", "mapproj", "maps", "sf", "raster", "spatstat", "dismo", "gstat"))

#Note: You need to have GDAL installed on your system
```

<h2>Reading the project files:</h2>
<p>

<p>

</p>
<p>
The <code>project_folder</code> needs to have the .rmd notebook and the .zip file inside of it. The script will unzip the file and recursively search all of the subfolders for .bil files:
</p>
```{r}
#Unzip the interpolation.zip file to a new directory called "interpolation" using the exdir parameter:
unzip(zipfile = "interpolation.zip", exdir = "interpolation")
```

<p>List the contents of the working directory in the console:</p>
```{r}
#Output the content of the interpolation folder to the console:
dir(path = "interpolation", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE)
```

<p>Using the recursive parameter of the dir function, we can list the contents of all the subfolders:</p>
```{r}
#Output the contents of the subfolders to the console:
dir(path = "interpolation/", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = TRUE)
```

<p>List only the .bil files:</p>
```{r}
#The pattern ="" parameter of the dir function allows us to list only file names matching a specified string pattern using regex
dir(path = "interpolation/", pattern ="\\.bil$", all.files = FALSE, full.names = FALSE, recursive = TRUE)
```

<p>Writing this list to a vector, we can store the names and locations of the .bil files in a string vector:</p>
```{r}
#Write the path of the files as a string vector to a new variable called "files" using "<-" R's assignment operator
files <- dir(path = "./interpolation/", pattern ="\\.bil$", all.files = FALSE, full.names = FALSE, recursive = TRUE)
```

<p>This the output is a vector list of the file locations:</p>
```{r}
#Calling the vector variable prints its contents in the console
files
```

<p>Concatenate <code>interpolation/</code> to the front of the files:</p>
```{r}
#paster() allows us concatenate a string to a vector of strings with one function call
paste("interpolation/", files, sep="", collapse = NULL)
```

<h2>Merging the raster data:</h2>
<p>Load the required <code>sf</code> and <code>raster</code> libraries:</p>
```{r message=FALSE, warning = FALSE}
library("sf")
library("raster")
```

<p>Use a for loop to add all of the rasters to an indexed list:</p>
```{r}
#We need to create an empty raster file
rasters <- NULL

#This for loop will run as many times as there are files:
for (val in 1:length(files)){
  rasters[[val]] <- raster(files[val])
}
```

<p>Merging the rasters together:
```{r}
#We will merger the first two rasters together:
merged_raster <- merge(rasters[[1]], rasters[[2]])

#An for loop will let us iterate through all of the raster files again
for (val in 3:length(files)){
  #We merge the current raster to the already merged rasters until we reach the end of the files[] 
  merged_raster <- merge(merged_raster, rasters[[val]])
}
```

<p>Print out the meta data for the merged_raster object:
```{r}
#Print() will output the default console output for a variable/object. In the case of a raster, some spatial statistics are displayed:
print(merged_raster)
```

<h2>Plotting the merged area:</h2>
<p>Load the required <code>ggplot</code> library (used for plotting), <code>ggplot</code> for color ramp generation, and set initial themes. The color palette in particular should be set to a palette suitable for representing elevation:</p>
```{r}
library("ggplot2")
library("RColorBrewer")

#Setting the default theme for ggplot2:
theme_set(theme_bw())
```

```{r}
#Setting a simple color theme with RColorBrewer:
palette <- rev(brewer.pal(n = 7, name = "RdYlGn"))

#Use the basic plot function to draw the merged raster on a map:
plot(merged_raster, main = "Plot of Area", xlab = "Longitude", ylab= "Latitude", col = palette, legend.args = list(text = 'Elevation (m)', side = 2))
```

<h2>Projecting the data:</h2>
<p>Create a string containing the name of the projection, ellipsoid type and a horizontal datum (NAD83 UTM Zone 11N):</p>
```{r}
#We can store the projection for the project as a string, as multiple functions may call it:
specified_projection <- "+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83"
```

<p><Code>raster</code> memory efficiency:
```{r}
#This Option sets the available memory fraction from the computer available for raster operations. A lower value will be more stable but offer slower performance (tested stable up to .3 in a machine with 16GB/3200Mhz RAM). Adjust this value as needed on a hardware by hardware basis:
rasterOptions(memfrac=.1)
```

<p>Projecting the data to (NAD83 UTM Zone 11N):</p>
```{r}
#This operation can take some time, depending on hardware:
projected_raster <- projectRaster(merged_raster, crs=specified_projection)
```

<p>Plot the projected data:</p>
```{r}
palette <- rev(brewer.pal(n = 7, name = "RdYlGn"))

plot(projected_raster, main = "Plot of Area (UTM Projected)", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation (m)', side = 2))
```

<p>Output the extent of the raster:</p>
```{r}
#extent() can be used to retrieve the extent of many types of variables, including rasters:
extent(projected_raster)
```

<h2>Creating the subject area polygon based on the center of the raster:</h2>
```{r}
#Set the side length to a variable equal to 100Km in metres (m), 1e+05
side_length <- 1e+05

#The x_origin and y_origin will be the center of the projected raster, we can find this by calling extent(projected_raster)[i], where i is the position of the xmin, ymax values... as shown in the previous console output from extent()):
x_origin <- ((extent(projected_raster)[1]+extent(projected_raster)[2])/2)+side_length/2
y_origin <- ((extent(projected_raster)[3]+extent(projected_raster)[4])/2)-side_length/2

#The x and y coords need to be a vector with four points, we then subtract the side length of the subject area to find the other corners"
x_coords <- c(x_origin, x_origin,  x_origin-side_length, x_origin-side_length)
y_coords <- c(y_origin, y_origin+side_length, y_origin+side_length, y_origin)

#The cbind() function will bind these to a new matrix:
subject_area <- cbind(x_coords, y_coords)

#Output the coordinates of the subject area polygon's vertexes:
subject_area
```
<p>Load the <code>sp</code> library and create the subject area polygon using the generated coordinates:</p>
```{r}
library(sp)

#We need to convert subject area coordinates to a polygon, we do so by passing them to the Polygon() function, then converting them to a SpatialPolygons sp object:
subject_area_polygon = Polygon(subject_area)
subject_area_polygons = Polygons(list(subject_area_polygon),1)
subject_area_SpatialPolygons = SpatialPolygons(list(subject_area_polygons))


plot(projected_raster, main = "Plot of Area", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation (m)', side = 2))

#We can plot subject_area_SpatialPolygons with the add = TRUE parameter to overlay it on top of the subject area raster:
plot(subject_area_SpatialPolygons, add = TRUE)
```

<h2>Plot the trimmed area:</h2>
<p>Trim the raster to the subject area polygon using the crop function:</p>
```{r}
#crop() will take an input raster as its first parameter and accept a SpatialPolygon as a clipping feature to crop the raster to
trimmed_subject_area <- crop(projected_raster, subject_area_SpatialPolygons)

plot(trimmed_subject_area, main = "Plot of Area", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation (m)', side = 2))

plot(subject_area_SpatialPolygons, add=TRUE)
```

```{r}
#hist() can be used to output a plot of the distribution of a raster:
hist(trimmed_subject_area, main="Histogram of Subject Area Elevation", xlab = "Elevation (m)", col= rev(brewer.pal(n = 5, name = "RdYlGn")), breaks = c(0,500,1000,1500,2000,2700))
```

<h2>Generate a field of random points:</h2>
```{r}
#This code block generates 1000 random points, with a minimum distance of 30m between points

#x and y are a our genesis points for the random field, generated as random deviates number between 0 and 1 to 15 decimal points 
x<- runif(1, 0, 1)
y<- runif(1, 0, 1)

#We create an empty matrix we can write the points to:
random_points <- cbind(x, y)

#This while loop is contingent on reaching 1000 valid sample points and will run until such values are found:
while(nrow(random_points) < 1000){
  #Each time the while loop runs, two new random x and y variables are generated for a new candidate random point:
  x_rand <- runif(1, 0, 1)
  y_rand <- runif(1, 0, 1)
  
  #We use a boolean variable to keep track of whether the minimum distance condition has been met by the candidate point:
  dist_check <- TRUE
  
  #A for loop will let us iterate through all of the valid points already found:
  for(j in 1:nrow(random_points)){
    #Simple coordinate inversion to calculate the distance between the current point in the list being checked, and the randomly generated candidate point:
    cur_dist <- sqrt(((random_points[j,1]-x_rand)^2)+((random_points[j,2]-y_rand)^2))
    
    #If the current distance is less than the minimum allowed distance, we set dist_check to false:
    if(cur_dist<0.003){
      dist_check <- FALSE
    }
  }
  
  #If the candidate point is far enough away from the rest of the points, dist_check will remain true, and we can bind it to our matrix of points:
  if(dist_check == TRUE){
    #rbind() allows you to append a new row to a matrix:
    random_points <- rbind(random_points,c(x_rand, y_rand))
  }
}

#A summary of some of the statistics of our randomly generated points:
summary(random_points)
```

<p>Scale and, transform and align the random points to the subject area:</p>
```{r}
#We need to scale the points to an appropriate level for the subject area
random_points_scaled <- random_points * 100000

#We Can recalculate an individual column easily with R, random_points_scaled[,1] corresponds to the x column
random_points_scaled[,1] <- random_points_scaled[,1] + x_origin  - 100000
random_points_scaled[,2] <- random_points_scaled[,2] + y_origin

#We convert these scaled points to a SpatialPoints object:
random_points_spatial <- SpatialPoints(coords = random_points_scaled)
random_points_geo <- SpatialPoints(random_points_scaled, proj4string=CRS(specified_projection)) 

plot(random_points_geo, main = "Plot of Random Points", xlab = "Eastings", ylab= "Northings",  pch=20)
plot(subject_area_SpatialPolygons, add=TRUE)

```

<p>Plot the trimmed subject area with the random points overlayed:</p>
```{r}
plot(trimmed_subject_area, main = "Plot of Area", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation (m)', side = 2))

plot(subject_area_SpatialPolygons, add=TRUE)
plot(random_points_geo, pch=20, cex=.5, add = TRUE)
```

<h2>Extract the raster value from the DEM at each random point:</h2>
```{r}
#extract() allows us to extract the value of a raster at the locations in SpatialPoints
extraction <- extract(trimmed_subject_area, random_points_geo, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, fun=NULL, na.rm=TRUE, layer, nl, df=FALSE, factors=FALSE)

#We can set a color ramp pallete with the maximum value of the extraction as the extent:
palette <- rev(colorRampPalette(c("red","orange","yellow", "yellowgreen", "greenyellow", "green", "darkgreen"))(max(extraction)))

plot(random_points_geo, main = "Plot of Area", xlab = "Eastings", ylab= "Northings", col = palette[cut(extraction, max(extraction))], pch=20) 
plot(subject_area_SpatialPolygons, add=TRUE)
```
<p>Print a histogram of the Random Point Elevations</p>
```{r}
hist(extraction, main = "Histogram of Random Point Elevation", col= rev(brewer.pal(n = 5, name = "RdYlGn")), breaks = c(0,500,1000,1500,2000,2700))
```

<h2> Creating Voronoi polygons:</h2>
<p>Load the <code>dismo</code> library:</p>
```{r message=FALSE, warning = FALSE}
library(dismo)
```

<p>Use the <code>voronoi</code> function to create a Voronoi tessellation:</p>
```{r}
vor <- voronoi(random_points_scaled)

palette <- rev(colorRampPalette(c("red","orange","yellow", "yellowgreen", "greenyellow", "green", "darkgreen"))(max(extraction)))

plot(vor, main = "Plot of Voronoi Polygons", xlab = "Eastings", ylab= "Northings", col = palette[cut(extraction, max(extraction))])
```
<p>Add the elevation data to a new vector in the Voronoi polygons data frame:</p>
```{r}
vor$elevation = extraction
```

<p>Obtain a summary of the Voronoi polygons attributes:</p>
```{r}
summary(vor)
```

<h2>Convert the Voronoi polygons to a raster data set:</h2>
```{r}
#Convert the vor_polygons to a SpatialPolygons data frame:
vor_polygons <- as(vor, "SpatialPolygonsDataFrame")
vor_polygons@data[,1] <- runif(nrow(vor_polygons))

#Create an empty raster to write the rasterization to and set its extent to the trimmed subject area:
vor_raster <- raster()
extent(vor_raster) <- extent(trimmed_subject_area)

#Perform the rasterization operation:
rasterization <- rasterize(vor_polygons, vor_raster, extraction, overlap='sum', res = res(trimmed_subject_area))

#Specify a projection for the new raster:
projection(rasterization) <- specified_projection

palette <- rev(colorRampPalette(c("red","orange","yellow", "yellowgreen", "greenyellow", "green", "darkgreen"))(max(extraction)))

#Output a simple plot of the results:
plot(rasterization,main = "Plot of Voronoi Polygons Rasterization", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation (m)', side = 2))
```
<p>Print a histogram of the Random Point Elevations:</p>
```{r}
hist(rasterization, main = "Histogram of Rasterized Voronoi Polygons from Random Points", xlab = "Elevation", col= rev(brewer.pal(n = 5, name = "RdYlGn")), breaks = c(0,500,1000,1500,2000,2700))
```
```{r}
summary(rasterization)
```

<h2>Calculating the difference between the source DEM and the Voronoi tessellation results:</h2>:
```{r}
#We need to resample the rasterization to be the same as the trimmed_subject area for the operations to work:
rasterization <- resample(rasterization, trimmed_subject_area)

#We can perform raster operations using primitive calls:
raster_difference <- rasterization - trimmed_subject_area 

#A summary of the difference between the two raster sets:
summary(raster_difference)
```

```{r}
#Simple plot of the difference:
palette <- rev(colorRampPalette(c("red","lightyellow", "blue"))(20))

plot(raster_difference, main = "Plot of Voronoi Polygons Raster Difference", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation Diff. (m)', side = 2))
```

<h2>Generating a inverse distance weighting interpolation (IDWI):</h2>
<p>Load the <code>spatstat</code> library:</p>
```{r message=FALSE, warning = FALSE}
library(spatstat)
```

<p>Generate the point pattern object:</p>
```{r}
#We use the ppp and idw functions from dismo to create the IDWI raster
ppp_random_points <- ppp(random_points_scaled[,1], random_points_scaled[,2], marks = extraction, window=owin(c(extent(random_points_scaled)[1], extent(random_points_scaled)[2]),c(extent(random_points_scaled)[3],extent(random_points_scaled)[4])))

random_points_scaled_idw <-idw(ppp_random_points)

random_points_scaled_idw_raster <- raster(random_points_scaled_idw)
```

<p>Plot the interpolated elevation map:</p>
```{r}
palette <- rev(colorRampPalette(c("red","orange","yellow", "yellowgreen", "greenyellow", "green", "darkgreen"))(max(extraction)))

plot(random_points_scaled_idw_raster, main = "Plot of Interpolated Points Raster", xlab = "Eastings", ylab= "Northings", col=palette, legend.args = list(text = 'Elevation (m)', side = 2))
```

```{r}
hist(random_points_scaled_idw_raster, main = "Histogram of IDWI from Random Points", xlab = "Elevation", col= rev(brewer.pal(n = 5, name = "RdYlGn")), breaks = c(0,500,1000,1500,2000,2500))
```

<h2>Calculating the difference between the subject area raster and the inverse distance weighting interpolation (IDWI) raster:</h2>:
```{r}
random_points_scaled_idw_raster <- resample(random_points_scaled_idw_raster, trimmed_subject_area)

raster_difference_idw <-random_points_scaled_idw_raster - trimmed_subject_area
```

```{r}
summary(raster_difference_idw)
```

```{r}
palette <- rev(colorRampPalette(c("red","lightyellow", "blue"))(20))

plot(raster_difference_idw, main = "Plot of IDWI Raster Difference", xlab = "Eastings", ylab= "Northings", col = palette, legend.args = list(text = 'Elevation Diff. (m)', side = 2))
```

<h2>Plotting maps side by side:</h2>

<p>We can add a second column to the parameters of our output to plot the two maps side by side:</p>
```{r}
#par() allows you to specify rows and columns for R's default output style
par(mfrow=c(1,2))

palette <- rev(colorRampPalette(c("red","orange","yellow", "yellowgreen", "greenyellow", "green", "darkgreen"))(max(extraction)))

plot(rasterization,main = "Voronoi Polygons", xlab = "Eastings", ylab= "Northings", col = palette, legend=FALSE, horizontal = TRUE)
plot(random_points_geo, pch=20, cex=.5, add = TRUE)

plot(random_points_scaled_idw_raster, main = "IDWI", col=palette, legend.args = list(text = 'Elevation (m)', side = 3), horizontal = TRUE)
plot(random_points_geo, pch=20, cex=.5, add = TRUE)
```

```{r}
par(mfrow=c(1,2))

palette <- rev(colorRampPalette(c("red","lightyellow", "blue"))(20))

plot(raster_difference, main = "Voronoi Polygons Difference",xlab = "Eastings", ylab= "Northings", col = palette, legend = FALSE, horizontal = TRUE)
plot(random_points_geo, pch=20, cex=.5, add = TRUE)

plot(raster_difference_idw, main = "IDWI Raster Difference", col = palette, legend.args = list(text = 'Elevation Diff. (m)', side = 3), horizontal = TRUE)
plot(random_points_geo, pch=20, cex=.5, add = TRUE)
```
<h2>Calculating the accuracy of the rasterized layers:</h2>
<i>"Note: you could quantify the total error in each interpolation by calculating the absolute value of each pixel in the error output (raster calculator), then summing the absolute pixel values (zonalstatistics as table with study area as zone). The largest sum is the most total error in the analysis, the smallest sum the least."</i> - Prof. C. Greene, <u><i>Introduction to Spatial Interpolation – TIN, IDW, NaturalNeighbours, Trend Surfaces, and Exploring Error</u></i>

<p id ="q1">Identify which interpolation method best approximated the continuous elevation surface and justify your answer.<p>
<p>We can calculate the sum of the absolute values of the raster cells using raster aggregate functions to calculate the total error, then compare:</p>
<p>For the Voronoi polygons rasterization:</p>
```{r}
#We apply the absolute function to the raster:
absolute_rasterized_voronoi_difference <- abs(raster_difference)
```

<p>Obtain the sum of the absolute values in the raster cells:</p>
```{r}
cellStats(absolute_rasterized_voronoi_difference, "sum")
```
<p id ="q1">Identify which interpolation method best approximated the continuous elevation surface and justify your answer.<p>
<p>For the IDWI rasterization:</p>
```{r}
#We apply the absolute function to the raster:
absolute_rasterized_idwi_difference <- abs(raster_difference_idw)
```

<p>Obtain the sum of the absolute values in the raster cells:</p>
```{r}
cellStats(absolute_rasterized_idwi_difference, "sum")
```
<p>At this particular subject area, and with this particular random sample rate and points, the Voronoi tessellation was the superior interpolation method- the calculated sum of the absolute difference for the Voronoi polygons is <code>1935755830(m)</code> < <code>2406133629(m)</code> the results for the interpolated distance weighting method.</p>

<p>In order to increase our certainty in this conclusion, it would be best to wrap the random sampling and Voronoi/IDWI rasterization and difference calculation in a loop, run in many times and then calculate the statistical mean error of the results.</p>


<h2>Add the layers to a leaflet map:</h2>
```{r warning=FALSE}
library(leaflet)

#In order to map the random points on the Leaflet map we will reproject them to WGS84:
random_points_geo <- spTransform(random_points_geo, CRS("+proj=longlat +datum=WGS84"))

#We specify the colors for our elevation mapping, using a conventional color scheme:
elevation_colors <- rev(c("white", "red", "orange", "yellow", "yellowgreen", "greenyellow", "green", "darkgreen"))
#seq() Allows you to generate a sequence based on a starting value, ending value and a width. For this particular map, a max elevation of 2700 will suffice, but this value could also be assigned using max(), for instance:
at_elevation <- seq(0, 2700, 200)
#colorBin allows you to create a bin based on these parameters: na.color = "#00000000" will ensure areas with no values will be transparent
cb_elevation <- colorBin(palette = elevation_colors , bins = at_elevation, domain = at_elevation, na.color = "#00000000")


#We specify the colors for our difference mapping, same procedure but different color names and manually assigned bin values:
difference_colors <- rev(c("red","lightyellow", "blue"))
at_difference <- c(-1500, -500, -100, 100, 500, 1500)
at_difference_display <- c(-1000, -500, -100, 100, 500, 1000)
cb_difference <- colorBin(palette = difference_colors , bins = at_difference, domain = at_difference, na.color = "#00000000")


#Create the leaflet object:
layeredMap <- leaflet() %>%
  #Add the source DEM to the map:
  addRasterImage(trimmed_subject_area, 
               colors = cb_elevation, 
               group = "Source DEM", 
               maxBytes= 9000000)  %>%
  
  #Add basemaps to the map:
  addProviderTiles(providers$OpenStreetMap.HOT,  group = "Open Street Map") %>%
  addProviderTiles(providers$OpenTopoMap,  group = "Open Topo Map") %>%
  addProviderTiles(providers$Esri.WorldImagery,  group = "Esri Imagery") %>%
  
  #Add the Voronoi results:
  addRasterImage(rasterization, 
                 colors = cb_elevation, #we use the Elevation color bins specified earlier  
                 group = "Voronoi Interpolation Results") %>% #A group name is required for modifying the legend
  
  #Add the IDWI results:
  addRasterImage(random_points_scaled_idw_raster, 
                 colors = cb_elevation, 
                 group = "IDWI Interpolation Results")%>%
  
  #Add the Voronoi Difference:
  addRasterImage(raster_difference, 
                 colors = cb_difference, #we use the Difference color bins specified earlier  
                 group = "Voronoi Difference", 
                 maxBytes= 9000000) %>%
  
  #Add the IDWI Difference:
  addRasterImage(raster_difference_idw, 
                 colors = cb_difference, 
                 group = "IDWI Difference", 
                 maxBytes= 9000000) %>%
  
  #Add the random sample point marker:
  addCircleMarkers(random_points_geo, 
                 lng = random_points_geo@coords[,1], #random_points_geo@coords[,1] will give us the lng column of our random points
                 lat= random_points_geo@coords[,2], #random_points_geo@coords[,2] will give us the lat column of our random points
                 radius = 1, 
                 weight = 1, 
                 color = "black", 
                 opacity = 1,
                 group = "Random Sample Points") %>% 
  
  #Add the elevation legend to the map:
  addLegend(pal = cb_elevation, 
            values = at_elevation, 
            title = "Elevation (m)") %>%
  
  #Add the elevation legend to the map:
  addLegend(pal = cb_difference, 
            values = at_difference_display, 
            title = "Difference (m)") %>%
  
  #Set the initial view for the map:
  setView(-120, 48, zoom = 8)%>%
  
  #Add the base layers (one must be visible to user):
  addLayersControl(
    baseGroups = c("Open Street Map","Open Topo Map", "Esri Imagery"),
    overlayGroups = c("Source DEM","Voronoi Interpolation Results","IDWI Interpolation Results", "Voronoi Difference","IDWI Difference", "Random Sample Points"),
    options = layersControlOptions(collapsed = FALSE, position = "bottomleft"))%>%
  
  #Hide all the layers:
  hideGroup(c("Voronoi Interpolation Results","IDWI Interpolation Results", "Voronoi Difference","IDWI Difference","Random Sample Points"))%>%
  
  #Force the Source DEM to be shown at initial load:
  showGroup("Source DEM")

#Output the map:
layeredMap
```
<script>
  document.getElementById('target').appendChild(document.getElementById(''))
</script>













