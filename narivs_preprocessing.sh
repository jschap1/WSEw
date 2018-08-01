# Command line commands for pre-processing narivs for cross sections analysis
# June 28, 2018

# ogrinfo -al -so  /Users/jschap/Documents/Data/HydroSHEDs/na_riv_15s/na_riv_15s.shp # get extent and projection

# Coordinate systems quick-reference:
# narivs, wgs84
# bathymetry pool, projected

# put in geographic coordinates
gdalwarp -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC/bath_pool_26/bath_1997_p26/w001001.adf" "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC/pool26_wgs.tif"

# get extent and projection
gdalinfo "/Users/jschap/Box Sync/Margulis_Research_Group/Jacob/UMBB/Data/UMESC/pool26_wgs.tif" -norat -noct

# Use extent to crop narivs
ogr2ogr -clipsrc -90.7511981 38.7913122 -90.1472941 39.0176724 /Users/jschap/Documents/Data/HydroSHEDs/na_riv_15s/Pools/nariv_cropped_p26.shp /Users/jschap/Documents/Data/HydroSHEDs/na_riv_15s/na_riv_15s.shp
