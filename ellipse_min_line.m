function [ x, y ] = ellipse_min_line( pos, theta, beta , ellipse )
%ELLIPSE_MIN_LIN Summary of this function goes here
%   Detailed explanation goes here

a = ellipse(1);
b = ellipse(2);

r = a*b/sqrt(b^2*cos(beta)*cos(beta)+a^2*sin(beta)*sin(beta));
R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];

vw_vec = r*[cos(beta); sin(beta)];

xy_vec = R*vw_vec;

x = [ pos(1), pos(1)+xy_vec(1)]; y=[pos(2), pos(2)+xy_vec(2)];

end

