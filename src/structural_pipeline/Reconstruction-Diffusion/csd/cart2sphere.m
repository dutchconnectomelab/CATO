function [r, theta, phi] = cart2sphere(x, y, z)
% [r, theta, phi] = CART2SPHERE(x, y, z) transforms cartesian coordinates 
% to spherical coordinates.
%
% given a point P with coordinates (x,y,z). r is the distance of P to the
% origin. theta is the angle between P and the Z-axis. phi is the angel
% between the projection of P onto the XY plane and the X-axis.

r = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./ r);
theta(r<=0) = 0;
phi = atan2(y, x);
