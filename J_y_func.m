function [ J_y ] = J_y_func( beta, theta, ellipse )
%J_X_FUNC Summary of this function goes here
%   Detailed explanation goes here

a = ellipse(1);
b = ellipse(2);

r = a*b./(sqrt(b^2*cos(beta).*cos(beta)+a^2*sin(beta).*sin(beta)));
l = r;

dbdt = - a*(1+(tan(pi-theta))^2)/(b*(1+cot(beta)^2));
dldb = a*b*(b^2-a^2)*(sin(beta)*cos(beta))*(b^2*(cos(beta))^2+a^2*(sin(beta))^2)^(-3/2);

J_y = -dbdt*(dldb*cos(beta-theta)-l*sin(beta-theta))-l*sin(beta-theta);

end

