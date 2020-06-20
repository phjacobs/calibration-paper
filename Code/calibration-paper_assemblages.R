## This code creates idealized SST and abundance relationships for
## planktic foraminifera, drawing upon previous work by Jonkers et al. 2019
## nature.com/articles/s41586-019-1230-3, including code by Jonkers
## accessible at https://zenodo.org/record/2638013.
## The final output of this code will be 160 synthetic assemblages
## collated into a single CSV as well separate CSVs per time interval.
## This code is available online at ZENODO LINK

## Section 0: Housekeeping
## 0.1 clear memory; set working directory

## 0.2. download Jonkers code and data:
## https://zenodo.org/record/2638013
## make_annual_fluxes_function.R
## make_annual_fluxes_assemblages.R
## shell flux data:  dat_sel.RDS
## standardize taxa: species_domains_compare.RDS

## 0.3. download ERSSTv5 SST data:
## https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.ersst.v5.html

## Section 1: Create Sed-Trap Assemblages
## use Jonkers et al. code to produce modern assemblages
#  source('make_annual_fluxes_function.R') sourced by next line
source('make_annual_fluxes_assemblages.R')
## should have TS.assem.RDS as output

## Section 2: Define SST preferences
## 2.1. extract SST for modern assemblages
## modified from Jonkers et al. 'extract_ERSSTv5.R' code
library(raster)
ERSST <- rotate(brick('sst.mnmean.nc'))

# get dates in ERSST
ERSST.date <- getZ(ERSST)
ERSST.month <- format(ERSST.date, '%Y-%m')
ERSST.year <- format(ERSST.date, '%Y')

# function to extract SST values using polygons
extractSST.fun <- function(dat, circle){
  tem <- extract(dat, circle, small = TRUE, weights= TRUE, df = TRUE)
  tem <- subset(tem, select = -ID)
  weights <- tem$weight
  vals <- subset(tem, select = -weight)
  if(nrow(tem)>1){
    apply(vals, 2, function(x) sum(x[!is.na(x)]*weights[!is.na(x)])/sum(weights[!is.na(x)]))}
  else{
    as.numeric(vals)
  }
} 
# historical SST sediment traps ####
# load polygons
trap_polys <- readRDS('trap_polys.RDS')
#trap.SST <- sapply(trap_polys, function(x) extractSST.fun(MSST, x))
#names(trap.SST) <- names(trap_polys)
#saveRDS(trap.SST, 'traps_ERSST_1854-1883.RDS')

# get mean temperature over collecting period of sediment trap ####
# needs annual assemblages (TS.assem) for length of the time series
get.mean.temp <- function(trap){
  begin <- format(as.Date(dat.sel[[which(names(dat.sel) == trap)]]$dat$open[1]), '%Y-%m')
  n.years <- nrow(TS.assem[[which(names(TS.assem) == trap)]])
  begin.indx <- which(ERSST.month == begin)
  end.indx <- (begin.indx + n.years*12)-1
  trap_poly <- trap_polys[[which(names(trap_polys) == trap)]]
  temperature <- extractSST.fun(ERSST[[begin.indx:end.indx]], trap_poly)
  mean.TS.SST <- mean(temperature, na.rm = TRUE)
  an.mean.TS.SST <- colMeans(matrix(temperature, nrow = 12), na.rm = TRUE)
  list(mean.TS.SST = mean.TS.SST, an.mean.TS.SST = an.mean.TS.SST)
}

dat.sel <- readRDS('dat_sel.RDS')
TS.assem <- readRDS('TS.assem.RDS')

TS.real.SST <- lapply(names(TS.assem), get.mean.temp)
names(TS.real.SST) <- names(TS.assem)
saveRDS(TS.real.SST, 'traps_ERSST_period.RDS')

## 2.2. create normal distribution for each taxa vs. SST


## Section 3: Generating Synthetic Samples
## 3.1. Virtual Distributions
## give VirtualSpecies Env raster and distribution

## 3.2. Virtual Samples
## restrict samples to a given ocean basin
## save Seed value to reuse for other intervals samples
## start with LGM (less ocean available than others)


## Section 4: Abundance-SST relationships
## 4.1. Linear Model of SST and Sed-TrapAbundance
## habitat suitability, taxa abundance, SSTs, other taxa's abundance


## 4.2. Predict Abundance for Samples


## Section 5: Save and Write Data to CSVs

#  all assemblages
#  modern assemblages
#  pre-industrial assemblages
#  mid-Holocene assemblages
#  Last Glacial Maximum assemblages