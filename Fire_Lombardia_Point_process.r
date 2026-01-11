### Use shape file and data Point Process
##
# File List 
# dataset: Incendi_2003.csv dataset - locations of fires in Lombardy in 2003
# WLR300.txt - defects detected on silicon wafers
# GaussBoaga Lombardy map: lombardia.shp, lombardia.dbf, lombardia.shx
# NUTS3 Austria map geographic coordinates: austrianuts3.shp, austrianuts3.dbf, austrianuts3.shx
#
# required libraries: maptools and spatstat (and dependencies)
# to check available packages, type library()
# for information on the individual library library(help = libraryname).
# Example: library(help = spatstat)
#
# PART 1: Reading shape files and displaying layers and points ------------
# Point pattern data
### Useful reference text for the sp and spdep libraries
### Bivand RS, Pebesma EJ, Gómez-Rubio V, 2008. Applied Spatial Data Analysis with R. Online book, Springer
#
#
require(maptools)   ### aggiunge il pacchetto maptools a sessione lavoro
require(rgdal)
# 
# reading shape file Lombardia in R (GaussBoaga)
lomb.poly <- readOGR("Lombardia.shp", verbose = TRUE) 
class(lomb.poly)
plot(lomb.poly)
slotNames(lomb.poly)
# lettura dati 
inc = incendi  <-  read.table("Incendi_2003.csv", header = T, sep = ";");   
str(incendi) 	      
# plot dati su mappa
points(incendi$EstN, incendi$Nord)
#
# PART 2: Example shape in geographic coordinates -------------------------
# area data
aus.poly <- readShapePoly("austrianuts3.shp", 
                        IDvar = "ID", verbose = TRUE) 			### read shape file in R
slotNames(aus.poly)
slot(aus.poly, "data")
aus.poly@data   # equivalent
a = slot(aus.poly, "polygons"); 
length(a); 
a[[1]]
plot(aus.poly, col = "green", border = "red")	### map region
centr = coordinates(aus.poly)
points(centr, col = "red")
text(centr[, 1], centr[, 2], aus.poly$ID, col = "blue")
xx <- aus.poly[aus.poly$ID == "AT13", ]    		# poligon Vienna
plot(xx, add = T, col = "blue")
#
# PART 3: PP: Using the shape to define a PPS window - package: spatstat ------------
require(spatstat) 
ppp <- incendi[, c("EstN", "Nord")]
ppp0 = as.ppp(ppp, W = lomb.poly);ppp0
# plot PP
plot(ppp0, cex = 0.5, main = , "Incendi in Lombardia nel 2003")
#
# PPS marked with extent of fire
ppp <- incendi[, c("EstN", "Nord", "Ettari")]
ppp0 = as.ppp(ppp, W = lomb.poly, mark = ppp$Ettari);ppp0
plot(ppp0, 
     main = , "Fires in Lombardy 2003 by extent"
)
#
d <- read.table("WLR300.txt", header = T) 	# read dataset in text format
str(d)						                  
#
W  <-  disc(96950, c(0, 0))			# defines a circular window for the process
pp <- ppp(d$X, d$Y, window = W); summary(pp)
plot(pp, cex = 0.6, main = "location events", cex.main  = 1.2, lwd = 3)

# PART 4: PP: Analysis of the CSR Hypothesis for Point Pattern Data - package: spatstat --------
# Test based on quadrat counts
qx <- quadratcount(pp, 4, 4)	### table 4x4
plot(qx)

te0 <- quadrat.test(pp, 4); te0					 ###test CSR Chi-squared
plot(te0, col = "red", cex = 1, lty = 1.5, lwd = 3, bg = gray(0.7)) ### counts, expected, residuals
plot(pp, pch = "+", cols = "green", cex = 1.5, lwd = 1.2, add = T)

qx10 <- quadratcount(pp, 10, 10)	### table 10x10
plot(qx10)
te10  <-  quadrat.test(pp, 10)	###test dispersion for 10x10
te10
# first evaluation clustering 
# histogram distances from nearest neighbor
hist(nndist(pp), main = "distance from nearest neighbor", xlab = "distance")
# test based on NN distribution
GGb <- Gest(pp, r = NULL, breaks = NULL, correction = "none");  
plot(GGb, raw ~ r, main = "CDF distance from NN", xlab = "r")
plot(GGb, theo~ r, add = T, col = 2, lty = 3)
legend(25000, 0.2,  legend = c("CDF empirical", "CDF teoric"), lty = c(1, 3), col = c(1, 2), bty = "n")
GGb	
#   test tipo G CDF corrected for border effect
GG <- Gest(pp, r = NULL, breaks = NULL);  
GG	
v <- plot(GG, cbind(rs, theo) ~ theo, main = "pp-plot")
#  Monte Carlo for test CSR
envpp <- envelope(pp, fun = Gest, nsim = 100, nrank = , verbose = TRUE, saveall = F)
win.graph()
aa <- plot(envpp, main = "MC", xlab = "y")
aa
# PART 5: PP: Kernel intensity estimation: intensity maps - package: spatstat --------
# intensity estimation
Z  <-  density.ppp(pp, 15000)
win.graph()
plot(Z, main = "kernel intensity"); 
plot(pp, add = T, cex = 0.6, col = "red")
persp(Z)
#
# Read shape file Lombardia
lomb.poly <- readShapePoly("lombardia.shp", verbose = TRUE) ### read shape file in R
# incendi_lombardia
ppp = incendi[, c("EstN", "Nord")]
ppp0 = as.ppp(ppp, W = lomb.poly); ppp0
plot(ppp0, cex = 0.5, main = , "Fires Lombardy in 2003")
# intensity estimation
Z  <-  density.ppp(ppp0, varcov = diag( c(var(ppp$EstN), var(ppp$Nord))/16)) 
win.graph()
plot(Z, main = "kernel intensity"); 
plot(ppp0, add = T, cex = 0.4)