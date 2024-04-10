#' Compute the arccosine of x
#'
#' This function computes the arccosine of the input value x.
#'
#' @param x Numeric vector.
#' @return Numeric vector containing the arccosine of each element of x.
#' @examples
#' Acos(0)
Acos = function(x)
{
  x = restrict(x, -1.0, +1.0)
  return(acos(x))
}

#' Compute the arcsine of x
#'
#' This function computes the arcsine of the input value x.
#'
#' @param x Numeric vector.
#' @return Numeric vector containing the arcsine of each element of x.
#' @examples
#' Asin(0)
#' @export
Asin = function(x)
{
  x = restrict(x, -1.0, +1.0)
  return(asin(x))
}

#' Compute the arc tangent of y/x
#'
#' This function computes the arc tangent of the input values y and x.
#'
#' @param y Numeric vector, the y-coordinate.
#' @param x Numeric vector, the x-coordinate.
#' @return Numeric vector containing the arctangent of y/x.
#' @export
#' @examples
#' Atan(0, 1)
Atan = function(y, x)
{
  # Special cases:
  if (x == 0.0)
  {
    if (0.0 < y)
      atan_yx = pi / 2.0
    else if (y < 0.0 )
      atan_yx = 3.0 * pi / 2.0
    else if ( y == 0.0 )
      atan_yx = 0.0
  }
  else if (y == 0.0)
  {
    if (0.0 < x)
      atan_yx = 0.0
    else if (x < 0.0)
      atan_yx = pi
  }
  else
  {

    abs_y = abs(y)
    abs_x = abs(x)
    theta = atan2(abs_y, abs_x)
    if (0.0 < x & 0.0 < y)
      atan_yx = theta
    else if (x < 0.0 & 0.0 < y)
      atan_yx = pi - theta
    else if (x < 0.0 & y < 0.0)
      atan_yx = pi + theta
    else if (0.0 < x & y < 0.0)
      atan_yx = 2.0 * pi - theta
  }
  return(atan_yx)
}

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

#' Normalize a vector
#'
#' This function normalizes the input vector v by dividing it by its L2 norm (Euclidean norm).
#'
#' @param v Numeric vector.
#' @return Numeric vector representing the normalized v.
#' @export
#' @examples
#' normalize(c(1,2,3))
normalize = function(v) {
  normalized_v = v / norm(v)
  return(normalized_v)
}

#' Normalize a matrix row-wise
#'
#' This function normalizes the rows of the input matrix x by applying the normalize function to each row.
#'
#' @param x Numeric matrix.
#' @return Numeric matrix with normalized rows.
#' @export
#' @examples
#' Normalize(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE))
Normalize = function(x)
{
  normalized_x = apply(x, 1, normalize)
  return(normalized_x)
}

#' Restrict a value to a specified range
#'
#' This function restricts the input value x to the specified range [lower, upper].
#'
#' @param x Numeric value to be restricted.
#' @param lower Numeric value representing the lower bound of the range.
#' @param upper Numeric value representing the upper bound of the range.
#' @return Numeric value restricted to the specified range.
#' @export
#' @examples
#' restrict(5, 0, 10)
restrict = function(x, lower, upper)
{
  x = max(x, lower)
  x = min(x, upper)
  return(x)
}

#' Cartesian to Spherical Coordinates Conversion
#'
#' Convert Cartesian coordinates to spherical coordinates.
#'
#' @param x A numeric vector of length 3 representing Cartesian coordinates (x, y, z).
#' @return A numeric vector of length 2 representing spherical coordinates (theta, phi).
#' @details
#' The Cartesian coordinates (x, y, z) are converted to spherical coordinates (theta, phi).
#' Theta represents the inclination angle (0 to pi), and phi represents the azimuth angle (0 to 2*pi).
#' @examples
#' cartesian_point <- matrix(c(1, 1, 1), nrow = 3)
#' Spherical_point <- cartesian_to_spherical(cartesian_point)
cartesian_to_spherical = function(x)
{
  theta = Acos(x[3])
  phi = Atan(x[2], x[1])
  spherical_coordinate = c(theta, phi)
  return(spherical_coordinate)
}

#' Convert Cartesian Coordinates to Spherical Coordinates
#'
#' This function converts Cartesian coordinates to spherical coordinates.
#'
#' @param x A matrix where each row represents a point in Cartesian coordinates.
#' @return A matrix where each row represents a point in spherical coordinates.
#' @details
#' The Cartesian coordinates (x, y, z) are converted to spherical coordinates (theta, phi).
#' Theta represents the inclination angle (0 to pi), and phi represents the azimuth angle (0 to 2*pi).
#' @examples
#' #example1
#' cartesian_points1 <- matrix(c(1, 1, 1,-1, 1, -1, 0, 0, 1),ncol = 3, byrow = TRUE)
#' Cartesian_to_Spherical(cartesian_points1)
#' #example2
#' cartesian_points2 <- matrix(c(4, 5, 9, -1, 0, 4, 5, 3, -1),ncol = 3, byrow = TRUE)
#' Cartesian_to_Spherical(cartesian_points2)
#' @export
Cartesian_to_Spherical = function(x)
{
  spherical_coordinates = t(apply(x, 1, cartesian_to_spherical))
  return(spherical_coordinates)
}

#' Compute the omega vector given spherical coordinates
#'
#' This function computes the omega vector given spherical coordinates theta and phi.
#'
#' @param theta Numeric value representing the polar angle.
#' @param phi Numeric value representing the azimuthal angle.
#' @return Numeric vector representing the omega vector.
#' @export
#' @examples
#' omega(pi/4, pi/3)
omega = function(theta, phi)
{
  omega_theta_phi = c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
  return(omega_theta_phi)
}

#' Convert Spherical Coordinates to Cartesian Coordinates
#'
#' This function converts spherical coordinates (theta, phi) to Cartesian coordinates.
#'
#' @param theta_phi Numeric vector of length 2 containing the spherical coordinates (theta, phi).
#' @return Numeric vector of length 3 containing the Cartesian coordinates (x, y, z).
spherical_to_cartesian = function(theta_phi)
{
  theta = theta_phi[1]
  phi = theta_phi[2]
  cartesian_coordinate = omega(theta, phi)
  return(cartesian_coordinate)
}

#' Convert Spherical Coordinates to Cartesian Coordinates for Multiple Points
#'
#' This function converts spherical coordinates (theta, phi) to Cartesian coordinates for multiple points.
#'
#' @param theta_phi Matrix where each row contains the spherical coordinates (theta, phi) of a point.
#' @return Matrix where each row contains the Cartesian coordinates (x, y, z) of a point.
#' @examples
#' theta_phi <- matrix(c(pi/4, pi/3, pi/6, pi/4), ncol = 2, byrow = TRUE)
#' Spherical_to_Cartesian(theta_phi)
#' @export
Spherical_to_Cartesian = function(theta_phi)
{
  cartesian_coordinates = t(apply(theta_phi, 1, spherical_to_cartesian))
  return(cartesian_coordinates)
}

#' Compute the cross product of two vectors
#'
#' This function computes the cross product of two input vectors u and v.
#'
#' @param u Numeric vector.
#' @param v Numeric vector.
#' @return Numeric vector representing the cross product of u and v.
#' @export
#' @examples
#' cross(c(1,0,0), c(0,1,0))
cross = function(u, v)
{
  cross_uv = c(u[2] * v[3] - u[3] * v[2],
               u[3] * v[1] - u[1] * v[3],
               u[1] * v[2] - u[2] * v[1])
  return(cross_uv)
}

#' Compute the normalized cross product of two vectors
#'
#' This function computes the cross product of two input vectors u and v,
#' normalizes the result, and returns the normalized vector.
#'
#' @param u Numeric vector.
#' @param v Numeric vector.
#' @return Numeric vector representing the normalized cross product of u and v.
#' @export
#' @examples
#' cross_normalized(c(1,0,0), c(0,1,0))
cross_normalized = function(u, v)
{
  cross_uv_normalized = normalize(cross(u, v))
  return(cross_uv_normalized)
}

#' Calculate Distance Between Two Vectors
#'
#' This function calculates the distance between two vectors using the dot product.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @return The distance between vectors x and y.
#' @examples
#' x <- c(1, 0, 0)
#' y <- c(0, 1, 0)
#' spherical_dist(x, y)
#' @export
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
#' @examples
#' edp(c(1,1,1))
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

#' Compute the exponential map on a Riemannian manifold
#'
#' This function computes the exponential map on a Riemannian manifold given a base point x and a tangent vector v.
#'
#' @param x Numeric vector representing the base point.
#' @param v Numeric vector representing the tangent vector.
#' @return Numeric vector representing the result of the exponential map.
#' @export
#' @examples
#' Exp(c(1,0,0), c(0,1,0))
Exp = function(x, v)
{
  if (sum(v^2) == 0)
    return(x)
  norm_v = norm2(v)
  Exp_x_v = cos(norm_v) * x + sin(norm_v) * v / norm_v
  return(Exp_x_v)
}

#' Compute points along the geodesic connecting two points on a Riemannian manifold
#'
#' This function computes points along the geodesic connecting two points p and q on a Riemannian manifold.
#'
#' @param t Time parameter for the geodesic path.
#' @param p Numeric vector representing the starting point.
#' @param q Numeric vector representing the ending point.
#' @param a Start time parameter.
#' @param b End time parameter.
#' @return Numeric vector representing a point along the geodesic path.
#' @export
#' @examples
#' geodesic(0.5, c(1,0,0), c(0,1,0), 0, 1)
geodesic = function(t, p, q, a, b)
{
  n = cross_normalized(p, q)
  w = cross_normalized(n, p)
  theta = spherical_dist(p, q) * (t - a) / (b - a)
  gamma = p * cos(theta) + w * sin(theta)
  return(gamma)
}

#' Compute points along the geodesic connecting two points on a Riemannian manifold at specified time points
#'
#' This function computes points along the geodesic connecting two points p and q on a Riemannian manifold at specified time points.
#'
#' @param t Numeric vector representing time points for the geodesic path.
#' @param p Numeric vector representing the starting point.
#' @param q Numeric vector representing the ending point.
#' @param a Start time parameter.
#' @param b End time parameter.
#' @return Numeric matrix representing points along the geodesic path at specified time points.
#' @export
#' @examples
#' Geodesic(c(0.25, 0.5, 0.75), c(1,0,0), c(0,1,0), 0, 1)
Geodesic = function(t, p, q, a, b)
{
  t = matrix(t, length(t), 1)
  gamma = t(apply(t, 1, geodesic, p, q, a, b))
  return(gamma)
}

#' Penalized Linear Spherical Spline
#'
#' This function fits a penalized linear spherical spline to the given data.
#'
#' This function estimates the optimal control points and knots for the given data, fits the model, and returns the optimized result.
#' Internally, gradient descent is used to minimize the loss of a given model, and a penalty term is added to control the complexity of the model.
#' Additionally, the BIC (Bayesian Information Criterion) value is calculated according to the model's complexity to provide information for model selection.
#' The function constructs piecewise curves using the piecewise_geodesic function and employs penalty terms to control the complexity of the model by updating control points and knots.
#' Its purpose is to find the optimal linear piecewise spline for the given data while controlling model complexity through penalty terms.
#' @param t A numeric vector representing the time or location.
#' @param y A matrix where each row represents a data point.
#' @param initial_control_points An optional matrix specifying initial control points. Default is NULL.
#' @param dimension An integer specifying the dimension of the spline.
#' @param initial_knots An optional numeric vector specifying initial knots. Default is NULL.
#' @param lambdas A numeric vector specifying the penalization parameters.
#' @param step_size A numeric value specifying the step size for optimization. Default is 0.01.
#' @param maxiter An integer specifying the maximum number of iterations. Default is 1000.
#' @param epsilon_iter A numeric value specifying the convergence criterion for iterations. Default is 1e-05.
#' @param jump_eps A numeric value specifying the threshold for pruning control points based on jump size. Default is 1e-02.
#' @param verbose A logical value indicating whether to print progress information. Default is FALSE.
#' @return A list containing the fit information.
#' @export
penalized_linear_spherical_spline= function(t, y, initial_control_points = NULL, dimension, initial_knots,
                                                  lambdas, step_size = 0.01, maxiter = 1000,
                                                  epsilon_iter = 1e-05, jump_eps = 1e-02, verbose = FALSE)
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
          control_point_tmp = Exp(control_points[j, ], -Rgrad_f * step)
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
        cat("jump size =", jump_size, "\n")
        cat("knotss =", as.numeric(knots), "\n")
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
#' # Define control points and knots
#' control_points <- matrix(c(2, 1, 0,   # Control point 1
#'                           3, 4, 2,   # Control point 2
#'                           5, 6, 3,   # Control point 3
#'                           7, 2, 1),  # Control point 4
#'                         nrow = 4, byrow = TRUE)
#' knots <- c(0, 1, 2, 3, 4)  # Knots indicating transitions
#' # Example of generating piecewise geodesic curve
#' t_example <- seq(0, 4, by = 0.1)
#' gamma_example <- piecewise_geodesic(t_example, control_points, knots)
#' # Plotting the piecewise geodesic curve
#' rgl.sphgrid(deggap = 90, col.long = "skyblue", col.lat = "skyblue")
#' spheres3d(x = 0, y = 0, z = 0, radius = 1, col = "grey",
#'          alpha = 0.05)
#' pch3d(rbind(v, w), col = "blue", cex = 0.3)
#' lines3d(gamma_example, col = "red", lty = 2)
piecewise_geodesic = function(t, control_points, knots)
{
  gamma = matrix(nrow = 0, ncol = 3)
  for (j in 1 : (nrow(control_points) - 1))
  {
    index_t = knots[j] <= t & t < knots[j + 1]
    piece_gamma = Geodesic(t[index_t], control_points[j, ], control_points[j + 1, ],
                           knots[j], knots[j + 1])
    gamma = rbind(gamma, piece_gamma)
  }
  return(gamma)
}

#' Calculate the Gradient of Loss Function for Linear Spline
#'
#' This function calculates the gradient of the loss function for a linear spline model.
#'
#' @param y Matrix of observed values.
#' @param t Vector of time or location values.
#' @param control_points Matrix of control points.
#' @param knots Vector of knots.
#' @param index Index of the control point.
#' @return Matrix containing the gradient of the loss function.
#' @export
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

#' Calculate Loss Function
#'
#' This function calculates the loss function based on the distance between observed values and predicted values.
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


#' Compute the Riemannian gradient penalty
#'
#' This function computes the Riemannian gradient penalty given control points, knots, and index.
#'
#' @param control_points Matrix representing the control points.
#' @param knots Numeric vector representing the knots.
#' @param index Index for computing the Riemannian gradient penalty.
#' @return A list containing Riemannian gradient and the exponential map of the negative Riemannian gradient.
#' @export
#' @examples
#' R_gradient_penalty(matrix(c(1,0,0,0,1,0,0,0,1), 3), c(0, 1, 2), 2)
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
  Exp_R_gradients = Exp(control_points[index, ], - R_gradients) / norm2(Exp(control_points[index, ], - R_gradients))
  return(list(R_grad = R_gradients, Exp_Rg = Exp_R_gradients))
}

#' Calculate Jump Vector for Linear Spline
#'
#' This function calculates the jump vector for a linear spline model.
#'
#' @param control_points Matrix of control points.
#' @param knots Vector of knots.
#' @return Jump vector.
#' @export
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
