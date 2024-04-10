# Define control points and knots
control_points <- matrix(c(2, 1, 0,   # Control point 1
                           3, 4, 2,   # Control point 2
                           5, 6, 3,   # Control point 3
                           7, 2, 1),  # Control point 4
                         nrow = 4, byrow = TRUE)
knots <- c(0, 1, 2, 3, 4)  # Knots indicating transitions

# Example of generating piecewise geodesic curve
t_example <- seq(0, 4, by = 0.1)
gamma_example <- piecewise_geodesic(t_example, control_points, knots)

# Plotting the piecewise geodesic curve
plot3D::scatter3D(gamma_example[, 1], gamma_example[, 2], gamma_example[, 3],
                  col = "blue", type = "l", main = "Example Piecewise Geodesic Curve",
                  xlab = "X", ylab = "Y", zlab = "Z")
