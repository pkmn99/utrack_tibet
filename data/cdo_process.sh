# reset lon to 0 to 360 and intert lat
cdo invertlat -sellonlatbox,0,360,-90,90 -setgrid,global_0.5 MODIS_landcover_05.nc processed/MODIS_landcover_05_clean.nc
cdo invertlat -sellonlatbox,0,360,-90,90 -setgrid,global_0.5 et_fraction.nc processed/et_fraction_clean.nc
