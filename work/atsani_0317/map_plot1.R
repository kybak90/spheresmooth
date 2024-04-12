library(rworldmap)
library(dplyr)
library(ggplot2)
library(geosphere)
#library(gpclib)
load("fit1_atsani.RData")
# World map
worldMap = getMap()
world.points = fortify(worldMap)
world.points$region = world.points$id

world.df = world.points[,c("long","lat","group", "region")]

colnames(world.df)[1] = "longitude"
colnames(world.df)[2] = "latitude"

points = data.frame(longitude = data_list$long,
                    latitude = data_list$lat,
                    time = spherical_data[,1],
                    central_press = data_list$central_press, 
                    wind_speed = data_list$wind_speed, 
                    grade = as.factor(data_list$grade), 
                    direction_50 = as.factor(data_list$direction_50),
                    long_radius_50 = data_list$long_radius_50,
                    short_radius_50 = data_list$short_radius_50,
                    direction_30 = as.factor(data_list$direction_30),
                    long_radius_30 = data_list$long_radius_30,
                    short_radius_30 = data_list$short_radius_30)




r3_best = Cartesian_to_Spherical(fit[[best_index]]$control_points)
r3_best
long_lat = r3_best*180 / pi

long_lat
points
long_lat_best = data.frame(latitude = 90-long_lat[, 1],
                           longitude = long_lat[,2])
long_lat_best

worldmap = ggplot() + 
   geom_polygon(data = world.df, aes(x = longitude, y = latitude, group = group),colour = "grey", fill = "antiquewhite") +
   scale_y_continuous(breaks = (-2:2) * 30) +
   scale_x_continuous(breaks = (-4:4) * 45) +
   geom_point(data = points, aes(x = longitude, y = latitude, color = central_press), cex = 0.8) +
   scale_color_gradient(low = "blue", high = "red") +
   geom_line(data = long_lat_best[1:3,], aes(x = longitude, y = latitude), col = 'black', lwd = 1) +
   geom_line(data = long_lat_best[3:5,], aes(x = longitude, y = latitude), col = 'black', lwd = 1) +
   geom_line(data = long_lat_best[5:6,], aes(x = longitude, y = latitude), col = 'black', lwd = 1) +
   geom_line(data = long_lat_best[6:8,], aes(x = longitude, y = latitude), col = 'black', lwd = 1) +
   geom_line(data = long_lat_best[8:9,], aes(x = longitude, y = latitude), col = 'black', lwd = 1) +
   coord_map("ortho", orientation=c(38, 120, 0))
worldmap


ggsave(worldmap, file='map_fitted.png', dpi = 300, width = 15, height = 15, units = "cm", bg = "transparent")



#####zoomed version
mar = 20
zoommap = ggplot() + 
   geom_polygon(data = world.df, aes(x = longitude, y = latitude, group = group),colour = "grey", fill = "antiquewhite") +
   scale_y_continuous(breaks = (-2:2) * 30) +
   scale_x_continuous(breaks = (-4:4) * 45) +
   geom_point(data = points, aes(x = longitude, y = latitude), cex = 0.8) +
   #scale_color_gradient(low = "blue", high = "red") +
   geom_line(data = long_lat_best[1:3,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
   geom_line(data = long_lat_best[3:5,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
   geom_line(data = long_lat_best[5:6,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
   geom_line(data = long_lat_best[6:7,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
   geom_line(data = long_lat_best[7:8,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
   geom_line(data = long_lat_best[8:9,], aes(x = longitude, y = latitude), col = 'red', lwd = 1) +
   coord_cartesian(xlim = c(min(long_lat_best$longitude) - mar, max(long_lat_best$longitude) + mar),
                   ylim = c(min(long_lat_best$latitude) - mar, max(long_lat_best$latitude) + mar)) +
   xlab('longitude') + ylab('latitude') +
   # for square-shape plot
   theme(axis.title=element_text(size = 7,face="bold"))

zoommap
ggsave(zoommap, file='map_fitted_atsani.png', dpi = 300, width = 17, height = 15, units = "cm", bg = "transparent")



#####zoomed version2
mar = 20
zoommap2 = ggplot() + 
   geom_polygon(data = world.df, aes(x = longitude, y = latitude, group = group),colour = "grey", fill = "antiquewhite") +
   scale_y_continuous(breaks = (-2:2) * 30) +
   scale_x_continuous(breaks = (-4:4) * 45) +
   geom_point(data = points, aes(x = longitude, y = latitude, color = central_press), cex = 0.8) +
   scale_color_gradient(low = "blue", high = "red") +
   geom_point(data = long_lat_best, aes(x = longitude, y = latitude), shape = 2, cex = 2.5) +
   coord_cartesian(xlim = c(min(long_lat_best$longitude) - mar, max(long_lat_best$longitude) + mar),
                   ylim = c(min(long_lat_best$latitude) - mar, max(long_lat_best$latitude) + mar)) +
   xlab('longitude') + ylab('latitude') +
   # for square-shape plot
   theme(axis.title=element_text(size = 7,face="bold"))


zoommap2

ggsave(zoommap2, file='map_fitted_atsani_central_press.png', dpi = 300, width = 17, height = 15, units = "cm", bg = "transparent")






# calculate value of jump size (penalty values)
control_points = fit[[best_index]]$control_points
knots = fit[[best_index]]$knots
number_penalty = nrow(control_points) - 2
jump_vector = matrix(0, number_penalty, 3)
velocity = vector(length = 8)
for (j in 1:8)
{
   theta1 = Acos(dot(control_points[j, ], control_points[j + 1, ]))
   velocity[j] = theta1/(knots[j+1] - knots[j])
}
velo = data.frame(knots = knots)
velo$v = 0
velo$v[1:8] = velocity

velocity_raw = vector(length = nrow(r3_data))
for (j in 1:60)
{
   theta1 = Acos(dot(r3_data[j, ], r3_data[j + 1, ]))
   velocity_raw[j] = theta1/(spherical_data[j+1,1] - spherical_data[j,1])
}



velo_raw  = data.frame(knots = spherical_data[, 1])
velo_raw$v = 0
velo_raw$v = velocity_raw

plot(velo_raw$knots, velo_raw$v, col = 'gray', xlab = "time point", ylab = "speed")
lines(velo$knots, velo$v)
a = rep(velo$knots[2:(nrow(velo)-1)], each = 2)
b = rep(velo$v[1:(nrow(velo)-1)], each = 2)

velo_step  = data.frame(knots = c(knots[1],a,knots[9]), v = c(b))
plot(velo_raw$knots, velo_raw$v, col = 'gray', xlab = "time point", ylab = "speed")
lines(velo_step$knots, velo_step$v)

############################################
plot(velo_raw$knots, velo_raw$v, col = 'gray', xlab = "time point", ylab = "speed")

lines(velo_step$knots[1:2], velo_step$v[1:2])
lines(velo_step$knots[2:3], velo_step$v[2:3], lty = 2)
lines(velo_step$knots[3:4], velo_step$v[3:4])
lines(velo_step$knots[4:5], velo_step$v[4:5], lty = 2)
lines(velo_step$knots[5:6], velo_step$v[5:6])
lines(velo_step$knots[6:7], velo_step$v[6:7], lty = 2)
lines(velo_step$knots[7:8], velo_step$v[7:8])
lines(velo_step$knots[8:9], velo_step$v[8:9], lty = 2)
lines(velo_step$knots[9:10], velo_step$v[9:10])
lines(velo_step$knots[10:11], velo_step$v[10:11], lty = 2)
lines(velo_step$knots[11:12], velo_step$v[11:12])
lines(velo_step$knots[12:13], velo_step$v[12:13], lty = 2)
lines(velo_step$knots[13:14], velo_step$v[13:14], col = 'blue')
lines(velo_step$knots[14:15], velo_step$v[14:15], lty = 2)
lines(velo_step$knots[15:16], velo_step$v[15:16], col = 'blue')
lines(c(velo_step$knots[16],knots[9]), c(velo_step$v[16],velo$v[9]), lty = 2)


points(velo_step$knots[c(1,3,5,7,9,11,13,15)], velo_step$v[c(1,3,5,7,9,11,13,15)], pch = 19,bg = 'black')
points(velo_step$knots[-c(1,3,5,7,9,11,13,15)], velo_step$v[-c(1,3,5,7,9,11,13,15)])
points(knots[9],velo$v[9], pch = 19, bg = 'black')

text(0.885, 0.5,"7", cex = 1.5)
text(0.96, 1.3, "8", cex = 1.5)
####################################################

points$year = round(data_list$time / 1000000)
points$day = round(data_list$time %% 1500000 / 100)
points$hour = data_list$time %% 100
points$date = points$year*1000 + points$day
test = as.Date(points$date, origin = "1972-04-27")
points$asdate = as.Date(test, "%Y-%m-%d")


# par(new=TRUE)
# plot(points$time, points$grade, type = 'l', axes = FALSE, xlab = "", ylab = "",
#      col = 'dark green')
# mtext("central press",side=4,col="dark green",line=4) 
# axis(4, col="dark green",col.axis="dark green",las=1)
#####zoomed version3
mar = 20
zoommap3 = ggplot() + 
   geom_polygon(data = world.df, aes(x = longitude, y = latitude, group = group),colour = "grey", fill = "antiquewhite") +
   scale_y_continuous(breaks = (-2:2) * 30) +
   scale_x_continuous(breaks = (-4:4) * 45) +
   geom_point(data = points, aes(x = longitude, y = latitude), cex = 0.8) +
   geom_segment(aes(x = long_lat_best$longitude[3], y = long_lat_best$latitude[3], xend = long_lat_jumps$longitude[1], 
                    yend = long_lat_jumps$latitude[1]), 
                arrow = arrow(),
                color='orange',size=0.03) +
   geom_point(data = long_lat_best, aes(x = longitude, y = latitude), shape = 2, cex = 2) +
   coord_cartesian(xlim = c(min(long_lat_best$longitude) - mar, max(long_lat_best$longitude) + mar),
                   ylim = c(min(long_lat_best$latitude) - mar, max(long_lat_best$latitude) + mar)) +
   xlab('longitude') + ylab('latitude') +
   # for square-shape plot
   theme(axis.title=element_text(size = 7,face="bold"))

zoommap3

ggsave(zoommap3, file='map_fitted_zoom1_central_press.png', dpi = 300, width = 15, height = 15, units = "cm", bg = "transparent")





plot(points$time, points$grade, type = 'l')
#axis(side = 1, at = fit[[best_index]]$knots)
abline(v = fit[[best_index]]$knots[1], col = 'red')
abline(v = fit[[best_index]]$knots[2], col = 'red')
abline(v = fit[[best_index]]$knots[3], col = 'red')
abline(v = fit[[best_index]]$knots[4], col = 'red')
abline(v = fit[[best_index]]$knots[5], col = 'red')
abline(v = fit[[best_index]]$knots[6], col = 'red')
abline(v = fit[[best_index]]$knots[7], col = 'red')


plot(points$asdate, points$wind_speed, type = 'l', xlab = 'time point', 
     ylab = 'grade')
#axis(side = 1, at = fit[[best_index]]$knots)
#abline(v = fit[[best_index]]$knots[1], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[2], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[3], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[4], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[5], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[6], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[7], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[8], col = 'red', lty = 2)


par(new=TRUE)
plot(points$time, points$central_press, type = 'l', axes = FALSE, xlab = "", ylab = "",
     col = 'dark green')
mtext("central press",side=4,col="dark green",line=4) 
axis(4, col="dark green",col.axis="dark green",las=1)

plot(points$time, points$long_radius_50, type = 'l',  xlab = 'time point', 
     ylab = 'longest radius for wind speed 50Kt')

points(points$time, points$long_radius_50, col = 'gray')
#points(points$time, points$wind_speed, col = ' gray')
#axis(side = 1, at = fit[[best_index]]$knots)
#abline(v = fit[[best_index]]$knots[1], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[2], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[3], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[4], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[5], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[6], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[7], col = 'red', lty = 2)



points$long_radius_50[is.na(points$long_radius_50)] = 0
plot(points$time, points$long_radius_50, type = 'l',  xlab = "time point", ylab = "longest radius for wind speed 50Kt")
#points(points$time, points$wind_speed, col = ' gray')
#axis(side = 1, at = fit[[best_index]]$knots)
#abline(v = fit[[best_index]]$knots[1], col = 'red', lty = 2)
abline(v = fit[[best_index]]$knots[2], col = 'red', lty = 2, lwd = 1.5)
abline(v = fit[[best_index]]$knots[3], col = 'red', lty = 2, lwd = 1.5)
abline(v = fit[[best_index]]$knots[4], col = 'blue', lty = 2, lwd = 1.5)
abline(v = fit[[best_index]]$knots[5], col = 'blue', lty = 2, lwd = 1.5)
abline(v = fit[[best_index]]$knots[6], col = 'blue', lty = 2, lwd = 1.5)
abline(v = fit[[best_index]]$knots[7], col = 'red', lty = 2, lwd = 1.5)
abline(v = fit[[best_index]]$knots[8], col = 'red', lty = 2, lwd = 1.5)

mtext("4", side = 3, at = c(0.63, 105), cex = 1.5)
mtext("5", side = 3, at = c(0.7, 105), cex = 1.5)
mtext("6", side = 3, at = c(0.77, 105), cex = 1.5)
