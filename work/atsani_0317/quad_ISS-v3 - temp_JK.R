####################################################################################
# Quadratic Intrinsic Spherical spline(QIS)
####################################################################################

# quadratic intrinsic spherical spline without penalty
fit_quadratic_spline = function(y, t, number_bezier = 3, initial_control_points = NULL,
                                step_size = 0.5, maxiter = 1000, 
                                epsilon_iter = 1e-05, verbose = FALSE)
{
   sample_size = length(t)
   J = number_bezier
   number_control_points = 2 * J + 1
   if (is.null(initial_control_points))
      control_points = y[seq(1, sample_size, length = number_control_points), ]
   else
      control_points = initial_control_points
   control_points = y[seq(1, sample_size, length = number_control_points), ]
   knots = s_knots_quantile(t, J + 1)
   gamma = quadratic_spline(t, control_points, knots)
   R = calculate_loss(y, gamma)
   temp_cp = control_points
   R_stored = Inf
   distance_R = rep(0, sample_size)
   for (iter in 1 : maxiter)
   {
      if (verbose)
         cat(iter, "th iteration runs \n")
      for (j in 1 : number_control_points)
      {
         cat(j)
         Rgrad_f = Rgradient_loss_quad_spline(y, t, control_points, knots, j)
         R_step = R
         step = step_size
         for (iter_step in 1 : 100)
         {  
            control_point_tmp = Exp(control_points[j, ], -Rgrad_f * step)
            temp_cp[j, ] = control_point_tmp / norm2(control_point_tmp)
            gamma = quadratic_spline(t, temp_cp, knots)
            R = calculate_loss(y, gamma)
            if (R_step >= R)
               break
            step = step / 2
         }
         control_points[j, ] = control_point_tmp
      }
      gamma = quadratic_spline(t, control_points, knots)
      R = calculate_loss(y, gamma)
      if (verbose)
         cat("R = ", R, "\n")
      if (abs(R - R_stored) < epsilon_iter)
         break
      R_stored = R
   }
   results = list(gamma = gamma, 
                  control_points = control_points)
   return(results)
}

# main code
# penalized linear sphrical spline
penalized_linear_spherical_spline_temp = function(t, y, initial_control_points = NULL, dimension, initial_knots,
                                                  lambdas, step_size = 0.01, maxiter = 1000,
                                                  epsilon_iter = 1e-05, jump_eps = 1e-02, verbose = FALSE)
{
   # browser()
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
   #if (sum(knots == t[31]) != 1)
   #{
   #   knots[min(which(knots > t[31])) -1] = t[31]
   #}
   Rlambda_stored = Inf
   distance_R = rep(0, sample_size)
   temp_cp = control_points
   # bic, dimension vector
   bic_list = rep(0, number_lambdas)
   #browser()
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
            #print(j)
            #if (j == 16)
               #brow
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
         # pruning step
         if (number_penalty > 0 & lambda > 0)
         {
            jump_size = rep(0, number_penalty)
            jumpC1 = jump_linear(control_points, knots)
            for (k in 1 : number_penalty)
            {
               weight = ((knots[k + 2] - knots[k + 1]) + (knots[k + 1] - knots[k])) / 2
               jump_size[k] = sum((weight * jumpC1[k, ])^2)
            }
            # for (k in 1 : number_penalty)
               # jump_size[k] = norm2(jump_vector[k, ])^2
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

# quadratic intrinsic spherical spline with only continuity condition
# the number of pieces of quadratic bezier curves are 
# the number of knots - 1, say 'J'
# the number of control points are 2J + 1
# each subinterval is [ knots[j-1], knots[1] )
quadratic_spline = function(t, control_points, knots)
{
   n = length(t)
   quad_iss = matrix(0, 0, 3)
   J = length(knots) - 1
   for (j in 1 : J)
   {
      I_j = knots[j] <= t & t < knots[j + 1]
      t_I_j = t[I_j]
      bezier_I_j = quadratic_bezier(t_I_j, 
                                    control_points[(2 * j - 1) : (2 * j + 1), ],
                                    knots[j : (j + 1)])
      quad_iss = rbind(quad_iss, bezier_I_j)
   }
   return(quad_iss)
}

# points for quadratic intrinsic spherical spline with continuity and C_1 conditions
# denote the number of pieces by 'J'
# the QISS is determined by (J + 2) control points
# the resulting number of points: 2J + 1
points_reflection = function(control_points, knots)
{
   J = length(knots) - 1
   # first, determine points 2J + 1 via reflections 
   if (J > 1)
   {
      points = matrix(0, 2 * J + 1, 3)
      points[1 : 3, ] = control_points[1 : 3, ]
      for (j in 1 : (J - 1))
      {
         v = points[2 * j, ]
         w = points[2 * j + 1, ]
         points[2 * j + 2, ] = reflection_knots(v, w, knots[j : (j + 2)])
         points[2 * j + 3, ] = control_points[j + 3, ]
      }
   }
   else
      points = control_points
   return(points)
}

# quadratic intrinsic spherical spline with continuity and C_1 conditions
# denote the number of pieces by 'J'
# the QISS is determined by (J + 2) control points
# each subinterval is [ knots[j-1], knots[1] )
quadratic_spline_reflection = function(t, control_points, knots)
{
   points = points_reflection(control_points, knots)
   quad_iss = quadratic_spline(t, points, knots)
   return(quad_iss)
}

# quadratic bezier curve
# three control points and two knots are required
# paramterization of time is (t - knots[j-1] / knots[j] - knots[j-1])
quadratic_bezier = function(t, control_points, knots)
{
   m = (t - knots[1]) / (knots[2] - knots[1])
   # range of time points [0, 1]
   p_t = Geodesic(m, control_points[1, ], control_points[2, ], 0, 1)
   q_t = Geodesic(m, control_points[2, ], control_points[3, ], 0, 1)
   c_t = calculate_c_t(p_t, q_t)
   R_1_c_t = calculate_R_s_c_t(c_t, 1 - m)
   R_c_t = calculate_R_s_c_t(c_t, m)
   quad_bezier = R_1_c_t * p_t + R_c_t * q_t
   return(quad_bezier)
}

# control_polygon
polygon_bezier = function(t, control_points)
{
   gamma = matrix(nrow = 0, ncol = 3)
   for (j in 1 : (nrow(control_points) - 1))
   {
      piece_gamma = Geodesic(t, control_points[j, ], control_points[j + 1, ], 0, 1)
      gamma = rbind(gamma, piece_gamma)
   }
   return(gamma)
}

# reflection of v in the space spanned by w at an real number a
# in fact, the a is determined by knots such that
# a = (knots[3] - knots[2]) / (knots[2] - knots[1]) dist(v, w)
reflection = function(v, w, a)
{
   b = dist(v, w)
   S_a_b = sin(a) / sin(b)
   reflection_v = (cos(a) + cos(b) * S_a_b) * w - S_a_b * v
   return(reflection_v)
}

# reflection of v in the space spanned by w at a
# the a is calculated by knots
reflection_knots = function(v, w, knots)
{
   b = dist(v, w)
   a = (knots[3] - knots[2]) / (knots[2] - knots[1]) * b
   reflection_v = reflection(v, w, a)
   return(reflection_v)
}

# calculate gradient of linear bezier curve
# two control points u, v determine the bezier curve
# there are two form of gradients for u, v
# by 'index' argument, it is determined
# output is returned in transpose form
gradient_linear_point = function(t, control_points, index)
{
   grad_linear = matrix(0, 3, 3)
   theta = dist(control_points[1, ], control_points[2, ])
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

# calculate gradient of quadratic bezier curve
# three control points u, v, w determine the bezier curve
# there are two three of gradients for u, v, w
# by 'index' argument, it is determined
gradient_quadratic_point = function(t, control_points, index)
{
   grad_quad = matrix(0, 3, 3)
   p = geodesic(t, control_points[1, ], control_points[2, ], 0, 1)
   q = geodesic(t, control_points[2, ], control_points[3, ], 0, 1)
   # for updating u, index = 1, so for p -> index, for q, index - 1
   # here when index - 1 = 0, not woriking
   # for updating v, index = 2, so for p -> index, for q, index - 1
   # for updating w, index = 3, so for p -> index, for q, index - 1
   # here when index = 3, not woriking; see 'gradient_linear_point'
   grad_p = gradient_linear_point(t, control_points[1 : 2, ], index)
   grad_q = gradient_linear_point(t, control_points[2 : 3, ], index - 1)
   c = dist(p, q)
   R_t = calculate_R_s(c, t)
   R_1_t = calculate_R_s(c, 1 - t)
   Q_t = calculate_Q_s(c, t)
   Q_1_t = calculate_Q_s(c, 1 - t)
   grad_R_t = Q_t * (grad_p %*% q + grad_q %*% p)
   grad_R_1_t = Q_1_t * (grad_p %*% q + grad_q %*% p)
   grad_quad = R_1_t * grad_p + R_t * grad_q + grad_R_1_t %*% t(p) + grad_R_t %*% t(q)
   return(grad_quad)
}

# for a point, 
# calculate Riemannian gradient of loss of linear bezier curve
Rgradient_loss_point_linear = function(y, t, control_points, index)
{
   grad_linear_bezier = gradient_linear_point(t, control_points, index)
   linear_bezier_t = Geodesic(t, control_points[1, ], control_points[2, ], 0, 1)
   phi = dist(y, linear_bezier_t)
   proj = calculate_projection_p(control_points[index, ], grad_linear_bezier %*% y)
   Aphi = calculate_Apsi(phi)
   Rgradient_linear_bezier = - Aphi * proj
   return(Rgradient_linear_bezier)
}

# for a point, 
# calculate Riemannian gradient of loss of quadratic bezier curve
Rgradient_loss_point_quad = function(y, t, control_points, index)
{
   grad_quad_bezier = gradient_quadratic_point(t, control_points, index)
   quad_bezier_t = quadratic_bezier(t, control_points, c(0, 1))
   phi = dist(y, quad_bezier_t)
   proj = calculate_projection_p(control_points[index, ], grad_quad_bezier %*% y)
   Aphi = calculate_Apsi(phi)
   Rgradient_quad_bezier = - Aphi * proj
   return(Rgradient_quad_bezier)
}

# for points, 
# calculate Riemannian gradient of loss of linear bezier curve
Rgradient_loss_linear = function(y, t, control_points, index)
{
   N = length(t)
   grad_linear_bezier = matrix(0, 1, 3)
   y = matrix(y, N, 3)
   for (n in 1 : N)
      grad_linear_bezier = grad_linear_bezier + t(Rgradient_loss_point_linear(y[n, ], t[n], control_points, index))
   return(grad_linear_bezier)
}

# for points, 
# calculate Riemannian gradient of loss of quadratic bezier curve
Rgradient_loss_quad = function(y, t, control_points, index)
{
   N = length(t)
   grad_quad_bezier = matrix(0, 1, 3)
   for (n in 1 : N)
      grad_quad_bezier = grad_quad_bezier + t(Rgradient_loss_point_quad(y[n, ], t[n], control_points, index))
   return(grad_quad_bezier)
}

# Riemannian gradient of loss for quadratic intrinsic spline
Rgradient_loss_quad_spline = function(y, t, control_points, knots, index)
{
   # divide subintervals via knots
   # first, check jth subinterval in 'index' 2j - 1 ~ 2j + 1
   K = nrow(control_points)
   J = length(knots) - 1
   ind_subintervals = seq(3, K, length = J)
   j = sum(!(index <= ind_subintervals)) + 1
   if (2 * j - 1 == index)
      sub_index = 1
   else if (2 * j == index)
      sub_index = 2
   else
      sub_index = 3
   sub_knots = knots[j : (j + 1)]
   sub_control_points = control_points[(2 * j - 1) : (2 * j + 1), ]
   I_j = knots[j] <= t & t < knots[j + 1]
   sub_y = y[I_j, ]
   sub_t = t[I_j]
   m_j = (sub_t - sub_knots[1]) / (sub_knots[2] - sub_knots[1])
   Rgrad = Rgradient_loss_quad(sub_y, m_j, sub_control_points, sub_index)
   if (sub_index == 3 & index != K)
   {
      j = j + 1
      sub_knots = knots[j : (j + 1)]
      sub_control_points = control_points[(2 * j - 1) : (2 * j + 1), ]
      I_j = knots[j] <= t & t < knots[j + 1]
      sub_y = y[I_j, ]
      sub_t = t[I_j]
      m_j = (sub_t - sub_knots[1]) / (sub_knots[2] - sub_knots[1])
      Rgrad = Rgrad + Rgradient_loss_quad(sub_y, m_j, sub_control_points, 1)
   }
   return(Rgrad)
}

# Riemannian gradient of loss for linear intrinsic spline
Rgradient_loss_linear_spline = function(y, t, control_points, knots, index)
{
   Rgrad = matrix(0, 1, 3)
   # divide subintervals via knots
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

# calculate the squared distance on the unit sphere
calculate_loss = function(y, gamma)
{
   n = nrow(y)
   distance = rep(0, n)
   for (i in 1 : n)
      distance[i] = dist(y[i, ], gamma[i, ])^2
   loss = 0.5 * sum(distance)
   return(loss)
}

# generate knots for quadratic intrinsic spherical spline
# here dimension is the number of knots
s_knots_quantile = function(x, dimension, tiny = 1e-5)
{
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

# calculate c(t) = dist(p(t), q(t))
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

# values at zero point of first derivative of linear bezier curve 
# in terms of time point
alpha_prime_zero = function(u, v)
{
   theta = dist(u, v)
   alpha_prime = theta * (- (cos(theta) / sin(theta)) * u + (1 / sin(theta)) * v)
   return(alpha_prime)
}

# values at one point of first derivative of linear bezier curve 
# in terms of time point
alpha_prime_one = function(u, v)
{
   theta = dist(u, v)
   alpha_prime = theta * (- (1 / sin(theta)) * u + (cos(theta) / sin(theta)) * v)
   return(alpha_prime)
}

# values at one point of first derivative of quadratic bezier curve 
# in terms of time point
beta_2prime_zero = function(u, v, w)
{
   a = dist(u, v)
   b = dist(v, w)
   tm = (alpha_prime_one(u, v) / sin(a)) - (alpha_prime_zero(u, v) / a)
   dev = (2 / sin(a)) * (a * alpha_prime_zero(v, w) - sin(a) * alpha_prime_zero(u, v)
                         + dot(alpha_prime_zero(v, w), u) * tm)
   return(dev)
}

# values at one point of first derivative of quadratic bezier curve 
# in terms of time point
beta_2prime_one = function(u, v, w)
{
   a = dist(u, v)
   b = dist(v, w)
   tm = (alpha_prime_zero(v, w) / sin(b)) - (alpha_prime_one(v, w) / b)
   dev = (2 / sin(b)) * (- b * alpha_prime_one(u, v) + sin(b) * alpha_prime_one(v, w)
                         + dot(alpha_prime_one(u, v), w) * tm)
   return(dev)
}

# calculate the Riemannian gradient of penalty function
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
         theta[k] = dist(control_points[index + k, ], control_points[index + k - 1, ])
      }
      # calculate the constant terms of penalty function
      a_j_minus = theta[1] / (sin(theta[1]) * delta[1])
      a_j = theta[2] / (sin(theta[2]) * delta[2])
      b_j = cos(theta[1]) * a_j
      b_j_minus = cos(theta[2]) * a_j_minus
      # difference of the corresponding penalty
      d = a_j_minus * control_points[index, ] - (b_j + b_j_minus) * control_points[index + 1, ] + a_j * control_points[index + 2, ]
      # values of derivetives
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
         theta[k] = dist(control_points[index + k - 1, ], control_points[index + k - 2 , ])
      }
      # calculate the constant terms of penalty function
      a_j_minus = theta[1] / (sin(theta[1]) * delta[1])
      a_j = theta[2] / (sin(theta[2]) * delta[2])
      b_j = cos(theta[1]) * a_j
      b_j_minus = cos(theta[2]) * a_j_minus
      # difference of the corresponding penalty
      d = a_j_minus * control_points[index - 1, ] - (b_j + b_j_minus) * control_points[index, ] + a_j * control_points[index + 1, ]
      # values of derivetives
      w1 = (cos(theta[2]) * theta[2] - sin(theta[2])) / sin(theta[2])^3
      w2 = (cos(theta[1]) * theta[1] - sin(theta[1])) / sin(theta[1])^3
      A1 = w1 * (- cos(theta[1]) * cos(theta[2]) * a_j_minus + dot(control_points[index - 1, ], control_points[index + 1, ]) * a_j_minus) - a_j
      B1 = theta[2] * cos(theta[2]) / sin(theta[2])
      A2 = w2 * (- cos(theta[2]) * cos(theta[1]) * a_j + dot(control_points[index - 1, ], control_points[index + 1, ]) * a_j) - a_j_minus
      B2 = theta[1] * cos(theta[1]) / sin(theta[1])
      # calculate the gradients
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
         theta[k] = dist(control_points[index + k - 2, ], control_points[index + k - 3, ])
      }
      # calculate the constant terms of penalty function
      a_j_minus = theta[1] / (sin(theta[1]) * delta[1])
      a_j = theta[2] / (sin(theta[2]) * delta[2])
      b_j = cos(theta[1]) * a_j
      b_j_minus = cos(theta[2]) * a_j_minus
      # difference of the corresponding penalty
      d = a_j_minus * control_points[index - 2, ] - (b_j + b_j_minus) * control_points[index - 1, ] + a_j * control_points[index, ]
      # values of derivetives
      A = (cos(theta[2]) * theta[2] - sin(theta[2])) / sin(theta[2])^3 * a_j_minus * dot(control_points[index - 2, ], control_points[index, ]) +
         (cos(theta[1]) * cos(theta[2]) * sin(theta[2]) - cos(theta[1]) * theta[2]) / sin(theta[2])^3 * a_j_minus - a_j
      B = theta[2] / sin(theta[2]) * a_j_minus
      # calculate the gradients
      grad = -(A * cos(theta[2]) + B * dot(control_points[index - 2, ], control_points[index, ])) * control_points[index, ] +
         A * control_points[index - 1, ] + B * control_points[index - 2, ]
      grad = grad / (norm2(d) * delta[2])
      R_gradients = R_gradients + grad
   }
   # calculate the Rimannian gradient and Exp(Rimannian gradient)
   R_gradients = t(R_gradients)
   Exp_R_gradients = Exp(control_points[index, ], - R_gradients) / norm2(Exp(control_points[index, ], - R_gradients))
   return(list(R_grad = R_gradients, Exp_Rg = Exp_R_gradients))
}

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

# calculate value of jump size (penalty values)
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

