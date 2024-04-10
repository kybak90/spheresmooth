cartesian_points1 <- matrix(c(1, 1, 1,
                             -1, 1, -1,
                             0, 0, 1),
                           ncol = 3, byrow = TRUE)

Cartesian_to_Spherical(cartesian_points1)


cartesian_points2 <- matrix(c(4, 5, 9,
                              -1, 0, 4,
                              5, 3, -1),
                            ncol = 3, byrow = TRUE)
Cartesian_to_Spherical(cartesian_points2)


