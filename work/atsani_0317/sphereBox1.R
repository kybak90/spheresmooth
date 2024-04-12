#  sphere box version 1.2

#  unrolling

#  transport
#  makes the parallel transport of xi along the geodesic 
#  gamma(t) = x cos t + w sin t, where x = gamma(0) and 
#  y = gamma(dist(x, y))
transport = function(xi, x, y)
{
   w = normalize(y - x * dot(x, y))
   d = dist(x, y)
   eta = xi + dot(w, xi) * (-sin(d) * x + (cos(d) - 1) * w)
   return(eta)
}

#  basic rotation operations

#  rotate_x
#  computes the rotation matrix with respect to x-axis
#  note that R stores columnwise
rotate_x = function(theta)
{
   rotate_x_theta = matrix(c(1, 0,           0,
                             0, cos(theta),  sin(theta),
                             0, -sin(theta), cos(theta)),
                           3, 3)
   return(rotate_x_theta)
}

#  rotate_y
#  computes the rotation matrix with respect to y-axis
#  this is equal to v_theta
#  note that R stores columnwise
rotate_y = function(theta)
{
   rotate_y_theta = matrix(c(cos(theta), 0, -sin(theta),
                             0,          1, 0,
                             sin(theta), 0, cos(theta)),
                           3, 3)
   return(rotate_y_theta)
}

#  rotate_z
#  computes the rotation matrix with respect to z-axis
#  this is equal to u_phi
#  note that R stores columnwise
rotate_z = function(phi)
{
   rotate_z_phi = matrix(c(cos(phi),  sin(phi), 0,
                           -sin(phi), cos(phi), 0,
                           0,         0,        1),
                         3, 3)
   return(rotate_z_phi)
}

#  variations of ratation matrices

#  rotate_a
#  computes the rotation mapping p and w to e3 and e1, respectively
rotate_a = function(p, q)
{
   r0 = rotate_p2e(p)
   w = normalize(q - p * dot(p, q))
   r0_w = r0 %*% w
   theta_phi = cartesian_to_spherical(r0_w)
   phi = theta_phi[2]
   r1 = rotate_z(-phi)
   ra = r1 %*% r0
   return(ra)
}

#  rotate_e2p
#  computes the rotation matrix which rotates e = (0, 0, 1) to p 
rotate_e2p = function(p)
{
   theta_phi = cartesian_to_spherical(p)
   theta = theta_phi[1]
   phi = theta_phi[2]
   r = rotate_z(phi) %*% rotate_y(theta)
   return(r)
}

#  rotate_geodesic
#  computes the rotation matrix mapping gamma(t) to e3, 
#  where gamma is the geodesic from p = gamma(a) to q = gamma(b)
rotate_geodesic = function(t, p, q, a, b, ra = NA)
{
   if (anyNA(ra))
      ra = rotate_a(p, q)
   rt = rotate_t(t, a, b, ra)
   return(rt)
}

#  rotate_p2e
#  computes the rotation matrix which rotates p to e = (0, 0, 1)
rotate_p2e = function(p)
{
   theta_phi = cartesian_to_spherical(p)
   theta = theta_phi[1]
   phi = theta_phi[2]
   r = rotate_y(-theta) %*% rotate_z(-phi)
   return(r)
}

#  rotate_p2q
#  computes the rotation matrix which rotates p to q
rotate_p2q = function(p, q)
{
   R_p2e = rotate_p2e(p)
   R_e2q = rotate_e2p(q)
   R_rotate_p2q = R_e2q %*% R_p2e
   return(R_rotate_p2q)
}

#  rotation_t
#  given ra, it computes the rotation matrix mapping gamma(t) to e3, 
#  where gamma is the geodesic from p = gamma(a) to q = gamma(b)
rotate_t = function(t, a, b, ra)
{
   theta = dist(p, q) * (t - a) / (b - a)
   rt = rotate_y(-theta) %*% ra
   return(rt)
}

#  basic operations on the sphere

#  cartesian_to_spherical
#  converts a cartesian coorinates x to the spherical coordinates (theta, phi)
cartesian_to_spherical = function(x)
{
   theta = Acos(x[3])
   phi = Atan(x[2], x[1])
   spherical_coordinate = c(theta, phi)
   return(spherical_coordinate)
}

#  Cartesian_to_Spherical
#  conerts a x-matrix to (theta, phi)-matrix
Cartesian_to_Spherical = function(x)
{
   spherical_coordinates = t(apply(x, 1, cartesian_to_spherical))
   return(spherical_coordinates)
}

#  cross
#  carries out the cross product of two vectors
cross = function(u, v)
{
   cross_uv = c(u[2] * v[3] - u[3] * v[2], 
                u[3] * v[1] - u[1] * v[3],
                u[1] * v[2] - u[2] * v[1])
   return(cross_uv)
}

#  cross_normalized
#  carries out the cross product of two vectors and normalize it
cross_normalized = function(u, v)
{
   cross_uv_normalized = normalize(cross(u, v))
   return(cross_uv_normalized)
}

#  dist
#  geodesic distance between two points on the sphere
dist = function(x, y)
{
   dist_x_y = Acos(dot(x, y))
   return(dist_x_y)
}

#  edp
#  projects a point on the sphere to a point on the plane 
#  preserving the distance from the northo pole
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

#  Exp
#  Exp maps v in Tx to a point in S3
Exp = function(x, v)
{
   if (sum(v^2) == 0)
      return(x)
   norm_v = norm2(v)
   Exp_x_v = cos(norm_v) * x + sin(norm_v) * v / norm_v
   return(Exp_x_v)
}

#  geodesic
#  computes gamma(t) which connects p and q with p = gamma(a) and q = gamma(b)
geodesic = function(t, p, q, a, b)
{
   n = cross_normalized(p, q)
   w = cross_normalized(n, p)
   theta = dist(p, q) * (t - a) / (b - a)
   gamma = p * cos(theta) + w * sin(theta)
   return(gamma)
}

#  Geodesic
#  computes gamma(t) which connects p and q with p = gamma(a) and q = gamma(b)
Geodesic = function(t, p, q, a, b)
{
   t = matrix(t, length(t), 1)
   gamma = t(apply(t, 1, geodesic, p, q, a, b))
   return(gamma)
}

#  Log
#  Log maps y in S3 to a vector in Tq
Log = function(x, y)
{
   if (norm2(x - y) < 1e-20)
      return(rep(0, 3))
   orthognal_complement = y - dot(x, y) * x
   orthognal_complement = orthognal_complement / norm2(orthognal_complement)
   Log_x_y = dist(x, y) * orthognal_complement
   return(Log_x_y)
}

#  omega
#  computes the euler angle representation of a point
omega = function(theta, phi)
{
   omega_theta_phi = c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
   return(omega_theta_phi)
}

#  spherical_to_cartesian
#  converts a (theta, phi) coorinates the x coordinates
spherical_to_cartesian = function(theta_phi)
{
   theta = theta_phi[1]
   phi = theta_phi[2]
   cartesian_coordinate = omega(theta, phi)
   return(cartesian_coordinate)
}

#  Cartesian_to_Spherical
#  conerts a (theta, phi)-matrix to the x-matrix
Spherical_to_Cartesian = function(theta_phi)
{
   cartesian_coordinates = t(apply(theta_phi, 1, spherical_to_cartesian))
   return(t(x))
}

#  weighted_average
#  computes the spherical weighted average sum wi pi
weighted_average = function(p, w,
                            max_iter = 1000,
                            epsilon = 1e-5)
{
   n = length(w)
   q = normalize(colSums(w * p))
   pstar = matrix(0, n, 3)
   rss = Inf
   for (iter in 1 : max_iter)
   {
      for (i in 1 : n)
         pstar[i, ] = Log(q, p[i, ])
      pstar_centroid = colSums(w * pstar)
      q_updated = Exp(q, pstar_centroid)
      if (norm2(q_updated - q) < epsilon)
         break
      q = q_updated
   }
   weighted_average = list(q = q, iter = iter)
   return(weighted_average)
}

# general functions

#  Acos
#  equals to acos when the argument is restriced to [-1, +1]
Acos = function(x)
{
   x = restrict(x, -1.0, +1.0)
   return(acos(x))
}

#  Asin
#  equals to asin when the argument is restriced to [-1, +1]
Asin = function(x)
{
   x = restrict(x, -1.0, +1.0)
   return(asin(x))
}

#  Atan
#  computes atan2(y, x) with correct angle
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
      #  We assume that atan2 is correct when both arguments are positive.
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

#  dot
#  dot product
dot = function(u, v)
{
   dot_u_v = sum(u * v)
   return(dot_u_v)
}

#  norm2
#  ell-2 norm
norm2 = function(u)
{
   ell2_norm_u = sqrt(dot(u, u))
   return(ell2_norm_u)
}

#  normalize
#  make a vector v to have norm one
normalize = function(v)
{
   normalized_v = v / norm2(v)
   return(normalized_v)
}

#  Normalize
#  make each row of x to have norm one
Normalize = function(x)
{
   normalized_x = apply(x, 1, normalize)
   return(normalized_x)
}

#  restrict
#  returns the value x restricted to [lower, upper]
restrict = function(x, lower, upper)
{
   x = max(x, lower)
   x = min(x, upper)
   return(x)
}

#  work notes
#  version 1.0 (2020.07.01)
#     sphereBox is launched
#  version 1.0 (2020.07.09)
#     several functions related to rotation are added
