This is most up-to-date version of the R codes for predicting unobserved bottom depth of river channels based on observations of their width and WSE. 

The workflow (1) for cross-section analysis is as follows:

find_centerlines
inputs: depths raster
outputs: river polyline

auto_transects
resample_polyline
auto_transects_helper_functions
inputs: depths raster, river polyline
outputs: transects data (distance, bed elevation, depth, smoothed, bankfull widths)

nice_cross_sections
fitting_sections_functions
fit_shape
inputs: transects data
outputs: shape parameters fit to "nice" cross sections

The idea is to extract many more transects than are actually used, since most of the transects are neither symmetric nor single-channeled.

***Workflow (2, current folder) uses the narivs centerlines, and takes a different approach to calculating width.

