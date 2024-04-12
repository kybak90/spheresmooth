library(rworldmap)
library(ggplot2)

load("./work/goni/goni_spherical.Rdata")

# spherical coordinates to Cartesian coordinates
goni_cartesian = Spherical_to_Cartesian(goni_spherical)

t = goni_spherical[, 1]
dimension = 15
initial_knots = knots_quantile(t, dimension = dimension)
lambda_seq = exp(seq(log(1e-07), log(1), length = 40))

goni_cartesian = Spherical_to_Cartesian(goni_spherical[, 2:3])

fit = penalized_linear_spherical_spline(t = t, y = goni_cartesian,
                                        dimension = dimension,
                                        initial_knots = initial_knots,
                                        lambdas = lambda_seq)
best_index = which.min(fit$bic_list)
best_index
# obtained control points for the piecewise geodesic curve
fit[[best_index]]$control_points

# world map visualization
worldMap = getMap()
world.points = fortify(worldMap)
world.points$region = world.points$id
world.df = world.points[, c("long","lat","group", "region")]
colnames(world.df)[1] = "longitude"
colnames(world.df)[2] = "latitude"

cp_best = Cartesian_to_Spherical(fit[[best_index]]$control_points)
cp_long_lat = cp_best * 180 / pi
cp_long_lat_df = data.frame(latitude = 90-cp_long_lat[, 1],
                           longitude = cp_long_lat[,2])

goni_spherical_df = data.frame(goni_spherical)
colnames(goni_spherical_df) = c("time", "latitude", "longitude")
goni_spherical_df$latitude = 90 - goni_spherical_df$latitude * 180 / pi
goni_spherical_df$longitude = goni_spherical_df$longitude * 180 / pi

worldmap = ggplot() +
  geom_polygon(data = world.df, aes(x = longitude, y = latitude, group = group),colour = "grey", fill = "antiquewhite") +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  geom_point(data = goni_spherical_df, aes(x = longitude, y = latitude), cex = 0.8) +
  geom_point(data = goni_spherical_df, aes(x = longitude, y = latitude), cex = 0.8) +
  geom_point(data = cp_long_lat_df, aes(x = longitude, y = latitude), shape = 23, col = 'blue', cex = 4) +
  geom_line(data = cp_long_lat_df[1:6,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
  geom_line(data = cp_long_lat_df[6:9,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
  geom_line(data = cp_long_lat_df[9:12,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
  coord_map("ortho", orientation=c(38, 120, 0))
worldmap
ggsave(worldmap, file='map_fitted_goni.png', dpi = 300, width = 15, height = 15, units = "cm", bg = "transparent")

# world map visualization (zoomed version)
mar = 20
zoommap = ggplot() +
  geom_polygon(data = world.df, aes(x = longitude, y = latitude, group = group),colour = "grey", fill = "antiquewhite") +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  geom_point(data = goni_spherical_df, aes(x = longitude, y = latitude), cex = 0.8) +
  geom_point(data = cp_long_lat_df, aes(x = longitude, y = latitude), shape = 23, col = 'blue', cex = 4) +
  geom_line(data = cp_long_lat_df[1:6,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
  geom_line(data = cp_long_lat_df[6:9,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
  geom_line(data = cp_long_lat_df[9:12,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
  coord_cartesian(xlim = c(min(cp_long_lat_df$longitude) - mar, max(cp_long_lat_df$longitude) + mar),
                  ylim = c(min(cp_long_lat_df$latitude) - mar, max(cp_long_lat_df$latitude) + mar)) +
  xlab('longitude') + ylab('latitude') +
  theme(axis.title=element_text(size = 7,face="bold"))
zoommap
ggsave(zoommap, file='map_fitted_goni_zoomed.png', dpi = 300, width = 17, height = 15, units = "cm", bg = "transparent")

