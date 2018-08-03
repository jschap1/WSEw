#' Resample Polyline
#'
#' Functions from agoldev for splitting polylines
#' Source: https://stackoverflow.com/questions/38700246/how-do-i-split-divide-polyline-shapefiles-into-equally-length-smaller-segments
#' The function splits a polyline into segments of a given length.
#' @param polyline: a spatial polyline data frame
#' @param interval_length: the length of the segments to split the lines into, in units of the polyline coordinates
#' @param add_original_points: whether or not the original points of the polyline should be added to the resulting line
#'                      if set FALSE, the resulting line will be shorter
#' @param add_final_point: whether or not the final point of the polyline should be added to the resulting line
#' @keywords resample, split, polyline
#' @export
#' @examples
#' resample_polyline(polyline, interval_length = section_length)

resample_polyline  = function(polyline, interval_length = 20, 
                              add_original_points = TRUE, add_final_point = FALSE) 
  {
  
  # transform input polyline
  linedf = data.frame(
    x  = polyline$x[1:nrow(polyline)-1],
    y  = polyline$y[1:nrow(polyline)-1],
    x2 = polyline$x[2:nrow(polyline)],
    y2 = polyline$y[2:nrow(polyline)]
  )
  
  # prepare output
  df = data.frame(
    x  = numeric(),
    y  = numeric()
  )
  
  residual_seg_length = 0
  for (i in 1:nrow(linedf)) {
    
    # for each line of the dataframe calculate segment length   
    v_seg      = linedf[i, ]
    seg_length = sqrt((v_seg$x - v_seg$x2) ^ 2 + (v_seg$y - v_seg$y2) ^ 2) 
    
    # create a vector of direction for the segment 
    v = c(v_seg$x2 - v_seg$x, v_seg$y2 - v_seg$y)
    
    # unit length
    u = c(v[1]  /  sqrt(v[1]  ^  2 + v[2]  ^  2), v[2]  /  sqrt(v[1]  ^  2 + v[2]  ^ 2))
    
    # calculate number of segment the segment is split into
    num_seg = floor((seg_length - residual_seg_length)  /  interval_length)
    
    # skip if next vertex is before interval_length
    if(num_seg >= 0) {      
      
      # add interpolated segments
      for (i in 0:(num_seg)) {
        df[nrow(df) + 1,] = c(
          v_seg$x  +  u[1] * residual_seg_length  +  u[1]  *  interval_length  *  i ,
          v_seg$y  +  u[2] * residual_seg_length  +  u[2]  *  interval_length  *  i
        )
      }
      
      # add original point (optional)
      if(add_original_points){
        df[nrow(df) + 1,] = c(
          v_seg$x2,
          v_seg$y2
        )
      }
      
    } else {
      
      # add original point (optional)
      if(add_original_points){
        df[nrow(df) + 1,] = c(
          v_seg$x2,
          v_seg$y2
        )
      }
      
      residual_seg_length = residual_seg_length - seg_length
      next()
      
    }
    
    # calculate residual segment length
    residual_seg_length = interval_length - ((seg_length - residual_seg_length) - (num_seg  *  interval_length))
    
  }
  
  # add final point (optional)
  if(add_final_point){
    df = rbind(df, data.frame(
      x = tail(polyline$x, n=1),
      y = tail(polyline$y, n=1)
    ))
  }
  
  return(df)
  
}
