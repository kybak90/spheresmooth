#' A spheresmooth package
#'
#' @description
#' Fitting a smooth path to a given set of noisy spherical data observed at known time points. It implements a piecewise geodesic curve fitting method on the unit sphere based on a velocity-based penalization scheme. The proposed approach is implemented using the Riemannian block coordinate descent algorithm. To understand the method and algorithm, one can refer to Bak, K. Y., Shin, J. K., & Koo, J. Y. (2023) <doi:10.1080/02664763.2022.2054962> for the case of order 1. Additionally, this package includes various functions necessary for handling spherical data.
#' @title Piecewise Geodesic Smoothing for Spherical Data
#' @name spheresmooth
NULL

#' Compute the dot product of two vectors
#'
#' This function computes the dot product of two input vectors u and v.
#'
#' @param u Numeric vector.
#' @param v Numeric vector.
#' @return Numeric value representing the dot product of u and v.
#' @export
#' @examples
#' dot(c(1,2,3), c(4,5,6))
dot = function(u, v) {
  dot_u_v = sum(u * v)
  return(dot_u_v)
}

#' Compute the L2 norm (Euclidean norm) of a vector
#'
#' This function computes the L2 norm (Euclidean norm) of the input vector u.
#'
#' @param u Numeric vector.
#' @return Numeric value representing the L2 norm of u.
#' @export
#' @examples
#' norm2(c(1,2,3))
norm2 = function(u)
{
  ell2_norm_u = sqrt(dot(u, u))
  return(ell2_norm_u)
}

# Normalize a vector
#'
#' This function normalizes the the input vector v by dividing its L2 norm (Euclidean norm).
#'
#' @param v Numeric vector.
#' @return Numeric vector with normalized.
#' @export
#' @examples
#' normalize_lower(1:6)
normalize_lower = function(v) {
  normalized_v = v / norm2(v)
  return(normalized_v)
}

#' Normalize a matrix row-wise
#'
#' This function normalizes the rows of the input matrix x by dividing each row by its L2 norm (Euclidean norm).
#'
#' @param x Numeric matrix.
#' @return Numeric matrix with normalized rows.
#' @export
#' @examples
#' normalize(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE))
normalize = function(x)
{
  normalized_x = apply(x, 1, normalize_lower)
  return(normalized_x)
}

#' Convert Cartesian coordinates to spherical coordinates
#'
#' This function converts Cartesian coordinates to spherical coordinates.
#'
#' @param x A matrix where each row represents a point in Cartesian coordinates.
#' @param byrow logical. If TRUE (the default) the matrix is filled by rows, otherwise the matrix is filled by columns.
#' @return A matrix where each row represents a point in spherical coordinates.
#' @details
#' The Cartesian coordinates (x, y, z) are converted to spherical coordinates (theta, phi).
#' Theta represents the inclination angle (0 to pi), and phi represents the azimuth angle (0 to 2*pi).
#' @examples
#' #example1
#' cartesian_points1 <- matrix(c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3),-1/sqrt(3), 1/sqrt(3), -1/sqrt(3)),
#'   ncol = 3, byrow = TRUE)
#' cartesian_to_spherical(cartesian_points1)
#' #example2
#' cartesian_points2 <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),ncol = 3, byrow = TRUE)
#' cartesian_to_spherical(cartesian_points2)
#' @export
cartesian_to_spherical = function(x, byrow = TRUE)
{
  if (is.vector(x))
    x = matrix(x, ncol = 3, byrow = byrow)
  theta = Acos(x[, 3])
  phi = Atan(x[, 2], x[, 1])
  spherical_coordinates = cbind(theta, phi)
  return(spherical_coordinates)
}



#' Convert spherical coordinates to Cartesian coordinates
#'
#' This function converts spherical coordinates (theta, phi) to Cartesian coordinates.
#'
#' @param theta_phi A matrix where each row contains the spherical coordinates (theta, phi) of a point.
#' @param byrow logical. If TRUE (the default) the matrix is filled by rows, otherwise the matrix is filled by columns.
#' @return A matrix where each row contains the Cartesian coordinates (x, y, z) of a point.
#' @examples
#' theta_phi <- matrix(c(pi/4, pi/3, pi/6, pi/4), ncol = 2, byrow = TRUE)
#' spherical_to_cartesian(theta_phi)
#' @export
spherical_to_cartesian = function(theta_phi, byrow = TRUE)
{
  if (is.vector(theta_phi))
    theta_phi = matrix(theta_phi, ncol = 2, byrow = byrow)
  cartesian_coordinates = omega(theta_phi[, 1], theta_phi[, 2])
  return(cartesian_coordinates)
}


#' Compute the cross product of two vectors
#'
#' This function computes the cross product of two input vectors u and v.
#'
#' @param u Numeric vector.
#' @param v Numeric vector.
#' @param normalize logical. If TRUE, returns the normalized vector of the cross product result.
#' @return Numeric vector representing the cross product of u and v.
#' @export
#' @examples
#' cross(c(1,0,0), c(0,1,0))
cross = function(u, v, normalize = FALSE)
{
  cross_uv = c(u[2] * v[3] - u[3] * v[2],
               u[3] * v[1] - u[1] * v[3],
               u[1] * v[2] - u[2] * v[1])
  if(normalize == TRUE)
    return(normalize_lower(cross_uv))
  return(cross_uv)
}

#' Calculate spherical distance between two vectors
#'
#' This function calculates the spherical distance between two vectors.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @return The distance between vectors x and y.
#' @export
#' @examples
#' x <- c(1, 0, 0)
#' y <- c(0, 1, 0)
#' spherical_dist(x, y)
spherical_dist = function(x, y)
{
  spherical_dist_x_y = Acos(dot(x, y))
  return(spherical_dist_x_y)
}

#' Compute the equal-distance projection of a point onto the xy plane
#'
#' This function computes the equal-distance projection of a point p onto the xy plane.
#'
#' @param p Numeric vector representing a point in Cartesian coordinates.
#' @return Numeric vector representing the equal-distance projection of p onto the xy plane.
#' @export
edp = function(p)
{
  theta_phi = cartesian_to_spherical(p)
  theta = theta_phi[1]
  phi = theta_phi[2]
  x = theta * cos(phi)
  y = theta * sin(phi)
  projection_p = c(x, y)
  return(projection_p)
}

#' Compute the exponential map on the unit sphere.
#'
#' This function computes the exponential map on the unit sphere given a base point x and a vector v.
#'
#' @param x Numeric vector representing the base point.
#' @param v Numeric vector representing a point.
#' @return Numeric vector representing the result of the exponential map.
#' @export
#' @examples
#' exp_map(c(0,0,1), c(1,1,0))
exp_map = function(x, v)
{
  if (sum(v^2) == 0)
    return(x)
  norm_v = norm2(v)
  Exp_x_v = cos(norm_v) * x + sin(norm_v) * v / norm_v
  return(Exp_x_v)
}

#' Compute the value of the geodesic curve connecting two points on the unit sphere for a given time point t
#'
#' This function computes points along the geodesic connecting two points p and q on the unit sphere.
#'
#' @param t Time parameter for the geodesic path.
#' @param p Numeric vector representing the starting point.
#' @param q Numeric vector representing the ending point.
#' @param a Start time parameter.
#' @param b End time parameter.
#' @return Numeric vector representing a point along the geodesic path.
#' @export
#' @examples
#' geodesic_lower(0.5, c(1,0,0), c(0,1,0), 0, 1)
geodesic_lower = function(t, p, q, a, b)
{
  n = cross(p, q, normalize = TRUE)
  w = cross(n, p, normalize = TRUE)
  theta = spherical_dist(p, q) * (t - a) / (b - a)
  gamma = p * cos(theta) + w * sin(theta)
  return(gamma)
}

#' Compute the value of the geodesic curve connecting two points on the unit sphere for a given set of time points t
#'
#' This function computes the value of the geodesic curve connecting two points p and q on the unit sphere at specified time points.
#'
#' @param t Numeric vector representing time points for the geodesic path.
#' @param p Numeric vector representing the starting point on the sphere.
#' @param q Numeric vector representing the ending point on the sphere.
#' @param a Start time parameter.
#' @param b End time parameter.
#' @return Numeric matrix representing points along the geodesic path at specified time points.
#' @export
#' @examples
#' geodesic(c(0.25, 0.5, 0.75), c(1,0,0), c(0,1,0), 0, 1)
geodesic = function(t, p, q, a, b)
{
  t = matrix(t, length(t), 1)
  gamma = t(apply(t, 1, geodesic_lower, p, q, a, b))
  return(gamma)
}

#' Generate knots for the piecewise geodesic curve based on the quantiles
#'
#' This generates a sequence of knots for a given set of time points based on the quantiles.
#'
#' @param x Numeric vector representing time points for the geodesic path.
#' @param dimension Numeric vector the number of knots.
#' @param tiny Numeric value representing a small constant that slightly expands the boundary.
#' @return Numeric vector representing knots sequence in the time domain.
#' @export
#' @examples
#' knots_quantile(seq(0, 1, length.out = 100), 10)
#' @importFrom stats quantile
knots_quantile = function(x, dimension, tiny = 1e-5)
{
  x = unique(x)
  dimension = max(dimension, 2)
  number_interior_knots = dimension - 2
  if (number_interior_knots > 0)
    probs = (1 : number_interior_knots) / (number_interior_knots + 1)
  else
    probs = NULL
  interior_knots = quantile(x[c(-1, -length(x))], probs, type = 3)
  knots = c(min(x) - tiny, interior_knots, max(x) + tiny)
  return(knots)
}

#' Penalized Linear Spherical Spline
#'
#' This function fits a penalized piecewise geodesic curve (linear spherical spline) to the given data.
#'
#' The goal is to find the optimal piecewise geodesic curve for the given spherical data while controlling model complexity through penalty terms.
#' This function computes the optimal control points and knots for the given data and returns the fitted result.
#' Internally, coordinate-wise gradient descent is used to minimize the loss function, and a penalty term is added to control the complexity of the model.
#' The BIC (Bayesian Information Criterion) value is calculated according to the model's complexity to provide information for model selection.
#' The function constructs piecewise curves using the piecewise_geodesic function and employs penalty terms to control the complexity of the model by updating control points and knots.
#' To see how to use the function in practical applications, refer to the README or https://github.com/kybak90/spheresmooth.
#'
#' @param t A numeric vector representing the time or location.
#' @param y A matrix where each row represents a data point on the sphere.
#' @param initial_control_points An optional matrix specifying initial control points. Default is NULL.
#' @param dimension An integer specifying the dimension of the spline.
#' @param initial_knots An optional numeric vector specifying initial knots. Default is NULL.
#' @param lambdas A numeric vector specifying the penalization parameters.
#' @param step_size A numeric value specifying the step size for optimization. Default is 1.
#' @param maxiter An integer specifying the maximum number of iterations. Default is 1000.
#' @param epsilon_iter A numeric value specifying the convergence criterion for iterations. Default is 1e-03.
#' @param jump_eps A numeric value specifying the threshold for pruning control points based on jump size. Default is 1e-04.
#' @param verbose A logical value indicating whether to print progress information. Default is FALSE.
#' @return A list containing the fitted result for each complexity parameter and BIC values for model selection. One might choose the element that corresponds to the minimum BIC values as illustrated in the example.
#' @export
#' @examples
#' library(sphereplot)
#' library(ggplot2)
#' library(sf)
#' library(rworldmap)
#' apw_cartesian = spherical_to_cartesian(apw_spherical[, 2:3])
#' t = apw_spherical[, 1]
#' dimension = 15
#' initial_knots = knots_quantile(t, dimension = dimension)
#' lambda_seq = exp(seq(log(1e-07), log(1), length = 40))
#' fit = penalized_linear_spherical_spline(t = t, y = apw_cartesian,
#'                                         dimension = dimension,
#'                                         initial_knots = initial_knots,
#'                                         lambdas = lambda_seq)
#' # choose a curve that minimizes the BIC
#' best_index = which.min(fit$bic_list)
#' best_index
#' # obtained control points for the piecewise geodesic curve
#' fit[[best_index]]$control_points
#'
#' worldMap = getMap()
#' worldMap_sf = st_as_sf(worldMap)
#' cp_best = cartesian_to_spherical(fit[[best_index]]$control_points)
#' cp_long_lat = cp_best * 180 / pi
#' cp_long_lat_df = data.frame(latitude = 90-cp_long_lat[, 1],
#'                             longitude = cp_long_lat[,2])
#' apw_spherical_df = data.frame(apw_spherical)
#' apw_spherical_df$latitude = 90 - apw_spherical_df$latitude * 180 / pi
#' apw_spherical_df$longitude = apw_spherical_df$longitude * 180 / pi
#' fitted_geodesic_curve = piecewise_geodesic(seq(0, 1, length = 2000),
#'                                            fit[[best_index]]$control_points,
#'                                            fit[[best_index]]$knots)
#'
#' fitted_cs = cartesian_to_spherical(fitted_geodesic_curve)
#' fitted_cs_long_lat = fitted_cs * 180 / pi
#' fitted_cs_long_lat_df = data.frame(latitude = 90 - fitted_cs_long_lat[, 1],
#'                                    longitude = fitted_cs_long_lat[, 2])
#'
#' apw_spherical_df_sf = st_as_sf(apw_spherical_df,
#'                                coords = c("longitude", "latitude"), crs = 4326)
#' cp_long_lat_df_sf = st_as_sf(cp_long_lat_df,
#'                              coords = c("longitude", "latitude"), crs = 4326)
#' fitted_cs_long_lat_df_sf = st_as_sf(fitted_cs_long_lat_df,
#'                                     coords = c("longitude", "latitude"), crs = 4326)
#'
#' # plot
#' worldmap = ggplot() +
#'   geom_sf(data = worldMap_sf, color = "grey", fill = "antiquewhite") +
#'   geom_sf(data = apw_spherical_df_sf, size = 0.8) +
#'   geom_sf(data = cp_long_lat_df_sf, color = "blue", shape = 23, size = 4) +
#'   geom_sf(data = fitted_cs_long_lat_df_sf, color = "red", size = 0.5) +
#'   xlab("longitude") +
#'   ylab("latitude") +
#'   scale_y_continuous(breaks = (-2:2) * 30) +
#'   scale_x_continuous(breaks = (-4:4) * 45) +
#'   coord_sf(crs = "+proj=ortho +lat_0=38 +lon_0=120 +y_0=0 +ellps=WGS84 +no_defs")
#' worldmap
penalized_linear_spherical_spline = function(t, y, initial_control_points = NULL, dimension, initial_knots,
                                                  lambdas, step_size = 1, maxiter = 1000,
                                                  epsilon_iter = 1e-03, jump_eps = 1e-04, verbose = FALSE)
{
  fit = list()
  sample_size = length(t)
  number_penalty = dimension - 2
  number_lambdas = length(lambdas)
  if (is.null(initial_control_points))
  {
    select_index = seq(1, sample_size, length = dimension)
    control_points = y[select_index, ]
  }
  else
    control_points = initial_control_points

  knots = initial_knots
  Rlambda_stored = Inf
  distance_R = rep(0, sample_size)
  temp_cp = control_points
  bic_list = rep(0, number_lambdas)
  dimension_list = rep(0, number_lambdas)
  gamma = piecewise_geodesic(t = t, control_points = temp_cp, knots = knots)
  Rlambda = calculate_loss(y, gamma)
  if (lambdas[1] > 0)
    Rlambda = Rlambda + lambdas[1] * sqrt(sum(rowSums(jump_linear(control_points, knots)^2)))
  for (lambda_index in 1 : number_lambdas)
  {
    lambda = lambdas[lambda_index]
    if (verbose)
      cat(lambda_index, "th lambda runs \n")
    for (iter in 1 : maxiter)
    {
      if (verbose)
        cat(iter, "th iteration runs \n")
      for (j in 1 : dimension)
      {
        Rgrad_loss = Rgradient_loss_linear_spline(y = y, t = t, control_points = control_points, knots = knots, index = j)
        R_pen = R_gradient_penalty(control_points, knots, j)$R_grad
        Rgrad_f = Rgrad_loss + lambda * R_pen
        R_step = Rlambda
        step = step_size
        for (iter_step in 1 : 100)
        {
          control_point_tmp = exp_map(control_points[j, ], -Rgrad_f * step)
          temp_cp[j, ] = control_point_tmp / norm2(control_point_tmp)
          gamma = piecewise_geodesic(t = t, control_points = temp_cp, knots = knots)
          Rlambda = calculate_loss(y, gamma)
          if (number_penalty > 0 & lambda > 0)
            Rlambda = Rlambda + lambda * sqrt(sum(rowSums(jump_linear(temp_cp, knots)^2)))
          if (R_step >= Rlambda)
            break
          step = step / 2
        }
        control_points[j, ] = temp_cp[j, ]
      }

      if (number_penalty > 0 & lambda > 0)
      {
        jump_size = rep(0, number_penalty)
        jumpC1 = jump_linear(control_points, knots)
        for (k in 1 : number_penalty)
        {
          weight = ((knots[k + 2] - knots[k + 1]) + (knots[k + 1] - knots[k])) / 2
          jump_size[k] = sum((weight * jumpC1[k, ])^2)
        }

        penalty_check = jump_size < jump_eps
        if (verbose)
        {
          cat("jump size =", jump_size, "\n")
          cat("knotss =", as.numeric(knots), "\n")
        }
        if (sum(penalty_check) > 0)
        {
          prune_index = which(penalty_check)
          control_points = control_points[-(prune_index + 1), ]
          knots = knots[-(prune_index + 1)]
          dimension = nrow(control_points)
          number_penalty = dimension - 2
          temp_cp = control_points
        }
      }
      gamma = piecewise_geodesic(t = t, control_points = temp_cp, knots = knots)
      Rlambda = calculate_loss(y, gamma)
      if (number_penalty > 0 & lambda > 0)
        Rlambda = Rlambda + lambda * sqrt(sum(rowSums(jump_linear(control_points, knots)^2)))
      if (verbose)
        cat("Rlambda =", Rlambda, "\n")
      if (abs(Rlambda - Rlambda_stored) < epsilon_iter)
        break
      Rlambda_stored = Rlambda
    }
    gamma = piecewise_geodesic(t = t, control_points = temp_cp, knots = knots)
    R = calculate_loss(y, gamma)
    dimension_list[lambda_index] = dimension
    bic_list[lambda_index] = sample_size * log(R) + 3 * log(sample_size) * dimension
    fit[[lambda_index]] = list(gamma = gamma,
                               control_points = control_points,
                               knots = knots,
                               dimension = dimension)
  }
  fit$dimension_list = dimension_list
  fit$bic_list = bic_list
  return(fit)
}


#' Piecewise Geodesic
#'
#' This function computes a piecewise geodesic path between control points.
#'
#' @param t A numeric vector representing the time or location.
#' @param control_points A matrix of control points where each row represents a control point.
#' @param knots A numeric vector of knot values.
#' @return A matrix containing the piecewise geodesic path.
#' @export
#' @details This function calculates the piecewise geodesic curve between control points based on the provided knots. The geodesic curve is computed segment by segment between adjacent control points. It interpolates the path between control points in a geodesic manner, ensuring the shortest path along the surface.
#' @examples
#' # `rgl` package and `sphereplot` pacakges are needed for the visualizaiton of the following example.
#' # Define control points and knots
#' library(rgl)
#' library(sphereplot)
#' control_points <- matrix(c(1, 0, 0,                         # Control point 1
#'                           1/sqrt(2), 1/sqrt(2), 0,          # Control point 2
#'                           -1/sqrt(3), 1/sqrt(3), 1/sqrt(3), # Control point 3
#'                           0, 0, 1),                         # Control point 4
#'                           nrow = 4, byrow = TRUE)
#' knots <- c(1, 2, 3, 3.5)  # Knots indicating transitions
#' # Example of generating piecewise geodesic curve
#' t_example <- seq(0, 4, by = 0.01)
#' gamma_example <- piecewise_geodesic(t_example, control_points, knots)
#' # Plotting the piecewise geodesic curve
#' rgl.sphgrid(deggap = 45, col.long = "skyblue", col.lat = "skyblue")
#' spheres3d(x = 0, y = 0, z = 0, radius = 1, col = "grey", alpha = 0.05)
#' pch3d(control_points, col = "blue", cex = 0.2, pch = 19)
#' lines3d(gamma_example, col = "red", lty = 1, lwd = 2)
piecewise_geodesic = function(t, control_points, knots)
{
  gamma = matrix(nrow = 0, ncol = 3)
  for (j in 1 : (nrow(control_points) - 1))
  {
    index_t = knots[j] <= t & t < knots[j + 1]
    piece_gamma = geodesic(t[index_t], control_points[j, ], control_points[j + 1, ],
                           knots[j], knots[j + 1])
    gamma = rbind(gamma, piece_gamma)
  }
  return(gamma)
}

#' Calculate Loss Function
#'
#' This function calculates the loss function based on the squared spherical distances between observed values and predicted values on the curve.
#'
#' @param y Matrix of observed values.
#' @param gamma Matrix of predicted values.
#' @return Loss value.
#' @export
calculate_loss = function(y, gamma)
{
  n = nrow(y)
  distance = rep(0, n)
  for (i in 1 : n)
    distance[i] = spherical_dist(y[i, ], gamma[i, ])^2
  loss = 0.5 * sum(distance)
  return(loss)
}

# Compute the arccosine of x
Acos = function(x)
{
  x = restrict(x, -1.0, +1.0)
  return(acos(x))
}

# Compute the arcsine of x
Asin = function(x)
{
  x = restrict(x, -1.0, +1.0)
  return(asin(x))
}

# Compute the arc tangent of y/x
Atan = function(y, x)
{
  # Handle special cases where x or y are zero
  atan_yx = numeric(length(y))

  zero_index = (x == 0)

  atan_yx[zero_index & y > 0] = pi / 2.0
  atan_yx[zero_index & y < 0] = 3.0 * pi / 2.0
  atan_yx[zero_index & y == 0] = 0.0

  nonzero_index = !zero_index & (y == 0)
  atan_yx[nonzero_index & x > 0] = 0.0
  atan_yx[nonzero_index & x < 0] = pi

  other_index = !(zero_index | y == 0)
  abs_y = abs(y[other_index])
  abs_x = abs(x[other_index])
  theta = atan2(abs_y, abs_x)

  atan_yx[other_index] = theta
  atan_yx[other_index & (x < 0 & y > 0)] = pi - theta[other_index & (x < 0 & y > 0)]
  atan_yx[other_index & (x < 0 & y < 0)] = pi + theta[other_index & (x < 0 & y < 0)]
  atan_yx[other_index & (x > 0 & y < 0)] = 2.0 * pi - theta[other_index & (x > 0 & y < 0)]

  return(atan_yx)
}

restrict = function(x, lower, upper)
{
  restrict_x = x
  restrict_x[x < lower] = lower
  restrict_x[x > upper] = upper
  return(restrict_x)
}

# Compute the omega vector given spherical coordinates
omega = function(theta, phi)
{
  omega_theta_phi = cbind(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
  return(omega_theta_phi)
}


# Compute the Riemannian gradient penalty
R_gradient_penalty = function(control_points, knots, index)
{
  # initial setting
  dimension = nrow(control_points)
  theta = rep(0, 2)
  delta = rep(0, 2)
  R_gradients = matrix(0, nrow = 3, ncol = 1)
  if (index < dimension - 1)
  {
    for (k in 1 : 2)
    {
      delta[k] = knots[index + k] - knots[index + k - 1]
      theta[k] = spherical_dist(control_points[index + k, ], control_points[index + k - 1, ])
    }

    a_j_minus = theta[1] / (sin(theta[1]) * delta[1])
    a_j = theta[2] / (sin(theta[2]) * delta[2])
    b_j = cos(theta[1]) * a_j
    b_j_minus = cos(theta[2]) * a_j_minus

    d = a_j_minus * control_points[index, ] - (b_j + b_j_minus) * control_points[index + 1, ] + a_j * control_points[index + 2, ]

    A = (cos(theta[1]) * theta[1] - sin(theta[1])) / sin(theta[1])^3 * a_j * dot(control_points[index, ], control_points[index + 2, ]) +
      (cos(theta[2]) * cos(theta[1]) * sin(theta[1]) - cos(theta[2]) * theta[1]) / sin(theta[1])^3 * a_j - a_j_minus
    B = theta[1] / sin(theta[1]) * a_j
    # calculate the gradients
    grad = -(A * cos(theta[1]) + B * dot(control_points[index, ], control_points[index + 2, ])) * control_points[index, ] +
      A * control_points[index + 1, ] + B * control_points[index + 2, ]
    grad = grad / (norm2(d) * delta[1])
    R_gradients = R_gradients + grad
  }
  if (index > 1 & index < dimension)
  {
    for (k in 1 : 2)
    {
      delta[k] = knots[index + k - 1] - knots[index + k - 2]
      theta[k] = spherical_dist(control_points[index + k - 1, ], control_points[index + k - 2 , ])
    }

    a_j_minus = theta[1] / (sin(theta[1]) * delta[1])
    a_j = theta[2] / (sin(theta[2]) * delta[2])
    b_j = cos(theta[1]) * a_j
    b_j_minus = cos(theta[2]) * a_j_minus

    d = a_j_minus * control_points[index - 1, ] - (b_j + b_j_minus) * control_points[index, ] + a_j * control_points[index + 1, ]

    w1 = (cos(theta[2]) * theta[2] - sin(theta[2])) / sin(theta[2])^3
    w2 = (cos(theta[1]) * theta[1] - sin(theta[1])) / sin(theta[1])^3
    A1 = w1 * (- cos(theta[1]) * cos(theta[2]) * a_j_minus + dot(control_points[index - 1, ], control_points[index + 1, ]) * a_j_minus) - a_j
    B1 = theta[2] * cos(theta[2]) / sin(theta[2])
    A2 = w2 * (- cos(theta[2]) * cos(theta[1]) * a_j + dot(control_points[index - 1, ], control_points[index + 1, ]) * a_j) - a_j_minus
    B2 = theta[1] * cos(theta[1]) / sin(theta[1])

    grad = (-(A1 * cos(theta[2]) + B1 * cos(theta[1])) * control_points[index, ] + A1 * control_points[index + 1, ] + B1 * control_points[index - 1, ]) / delta[2] +
      (-(A2 * cos(theta[1]) + B2 * cos(theta[2])) * control_points[index, ] + A2 * control_points[index - 1, ] + B2 * control_points[index + 1, ]) / delta[1]
    grad = grad / (norm2(d))
    R_gradients = R_gradients + grad
  }
  if (index > 2)
  {
    for (k in 1 : 2)
    {
      delta[k] = knots[index + k - 2] - knots[index + k - 3]
      theta[k] = spherical_dist(control_points[index + k - 2, ], control_points[index + k - 3, ])
    }

    a_j_minus = theta[1] / (sin(theta[1]) * delta[1])
    a_j = theta[2] / (sin(theta[2]) * delta[2])
    b_j = cos(theta[1]) * a_j
    b_j_minus = cos(theta[2]) * a_j_minus

    d = a_j_minus * control_points[index - 2, ] - (b_j + b_j_minus) * control_points[index - 1, ] + a_j * control_points[index, ]

    A = (cos(theta[2]) * theta[2] - sin(theta[2])) / sin(theta[2])^3 * a_j_minus * dot(control_points[index - 2, ], control_points[index, ]) +
      (cos(theta[1]) * cos(theta[2]) * sin(theta[2]) - cos(theta[1]) * theta[2]) / sin(theta[2])^3 * a_j_minus - a_j
    B = theta[2] / sin(theta[2]) * a_j_minus

    grad = -(A * cos(theta[2]) + B * dot(control_points[index - 2, ], control_points[index, ])) * control_points[index, ] +
      A * control_points[index - 1, ] + B * control_points[index - 2, ]
    grad = grad / (norm2(d) * delta[2])
    R_gradients = R_gradients + grad
  }

  R_gradients = t(R_gradients)
  Exp_R_gradients = exp_map(control_points[index, ], - R_gradients) / norm2(exp_map(control_points[index, ], - R_gradients))
  return(list(R_grad = R_gradients, Exp_Rg = Exp_R_gradients))
}

# Calculate Jump Vector for Linear Spline
jump_linear = function(control_points, knots)
{
  number_penalty = nrow(control_points) - 2
  jump_vector = matrix(0, number_penalty, 3)
  for (j in 1 : number_penalty)
  {
    theta1 = Acos(dot(control_points[j, ], control_points[j + 1, ]))
    theta2 = Acos(dot(control_points[j + 1, ], control_points[j + 2, ]))
    a = theta2 / sin(theta2) / (knots[j + 2] - knots[j + 1])
    c = theta1 / sin(theta1) / (knots[j + 1] - knots[j])
    b1 = cos(theta2) * a
    b2 = cos(theta1) * c
    jump_vector[j, ] = a * control_points[j + 2, ] -
      (b1 + b2) * control_points[j + 1, ] + c * control_points[j, ]
  }
  return(jump_vector)
}

# Calculate the Gradient of Loss Function for Linear Spline Curve
Rgradient_loss_linear_spline = function(y, t, control_points, knots, index)
{
  Rgrad = matrix(0, 1, 3)
  J = nrow(control_points)
  if (index > 1)
  {
    sub_knots = knots[(index - 1) : index]
    sub_control_points = control_points[(index - 1) : index, ]
    I_j = knots[index - 1] <= t & t < knots[index]
    sub_y = y[I_j, ]
    sub_t = t[I_j]
    m_j = (sub_t - sub_knots[1]) / (sub_knots[2] - sub_knots[1])
    Rgrad = Rgradient_loss_linear(sub_y, m_j, sub_control_points, 2)
  }
  if (index < J)
  {
    sub_knots = knots[index : (index + 1)]
    sub_control_points = control_points[index : (index + 1), ]
    I_j = knots[index] <= t & t < knots[index + 1]
    sub_y = y[I_j, ]
    sub_t = t[I_j]
    m_j = (sub_t - sub_knots[1]) / (sub_knots[2] - sub_knots[1])
    Rgrad = Rgrad + Rgradient_loss_linear(sub_y, m_j, sub_control_points, 1)
  }
  return(Rgrad)
}

gradient_linear_point = function(t, control_points, index)
{
  grad_linear = matrix(0, 3, 3)
  theta = spherical_dist(control_points[1, ], control_points[2, ])
  Q_t = calculate_Q_s(theta, t)
  Q_1_t = calculate_Q_s(theta, 1 - t)
  if (index == 1)
  {
    R_1_t = calculate_R_s(theta, 1 - t)
    grad_linear = R_1_t * diag(1, 3) + Q_1_t * outer(control_points[2, ], control_points[1, ]) +
      Q_t * outer(control_points[2, ], control_points[2, ])
  }
  else if (index == 2)
  {
    R_t = calculate_R_s(theta, t)
    grad_linear = R_t * diag(1, 3) + Q_1_t * outer(control_points[1, ], control_points[1, ]) +
      Q_t * outer(control_points[1, ], control_points[2, ])
  }
  else
    grad_linear = matrix(0, 3, 3)
  return(grad_linear)
}

Rgradient_loss_point_linear = function(y, t, control_points, index)
{
  grad_linear_bezier = gradient_linear_point(t, control_points, index)
  linear_bezier_t = geodesic(t, control_points[1, ], control_points[2, ], 0, 1)
  phi = spherical_dist(y, linear_bezier_t)
  proj = calculate_projection_p(control_points[index, ], grad_linear_bezier %*% y)
  Aphi = calculate_Apsi(phi)
  Rgradient_linear_bezier = - Aphi * proj
  return(Rgradient_linear_bezier)
}

Rgradient_loss_linear = function(y, t, control_points, index)
{
  N = length(t)
  grad_linear_bezier = matrix(0, 1, 3)
  y = matrix(y, N, 3)
  for (n in 1 : N)
    grad_linear_bezier = grad_linear_bezier + t(Rgradient_loss_point_linear(y[n, ], t[n], control_points, index))
  return(grad_linear_bezier)
}

# calculate c(t)
#' @importFrom stats dist
calculate_c_t = function(p_t, q_t)
{
  c_t = rep(0, nrow(p_t))
  for (i in 1 : nrow(p_t))
    c_t[i] = dist(p_t[i, ], q_t[i, ])
  return(c_t)
}

# calculate R_s(c(t))
calculate_R_s_c_t = function(c_t, s)
{
  R_s_c_t = rep(0, length(s))
  for (i in 1 : length(s))
    R_s_c_t[i] = calculate_R_s(c_t[i], s[i])
  return(R_s_c_t)
}

# R(theta, s) = sin(s theta) / sin(theta)
calculate_R_s = function(theta, s)
{
  if (theta == 0)
    return(0)
  else
    return(sin(theta * s) / sin(theta))
}

# calculate Q_s(theta)
calculate_Q_s = function(theta, s)
{
  value_Qs = (sin(s * theta) * cos(theta) - s * cos(s * theta) * sin(theta)) / sin(theta)^3
  return(value_Qs)
}

# psi / sin(psi)
calculate_Apsi = function(psi)
{
  if (sum(psi^2) == 0)
    return(0)
  else
    return(psi / sin(psi))
}

# projection y on sphere onto tangent plane at p
calculate_projection_p = function(p, y)
{
  proj_y = (y - p * dot(p, y)) #/ norm2(y - p * dot(p, y))
  return(proj_y)
}
