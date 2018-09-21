setwd("/Users/jschap/Desktop/Cross_Sections")

load("Data/Processed/processed_xs_data.rda")

plot(WSE~w, rWSEw[[3000]])
plot(WSE~w, xWSEw[[3000]])
plot(cross_sections$x[3000][[1]], cross_sections$b[3000][[1]], type="l")
