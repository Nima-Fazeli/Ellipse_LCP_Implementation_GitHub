function [ Y_min, beta_c ] = min_point_ellipse( pos, ori, ellipse, res )
%MIN_POINT_ELLIPSE Summary of this function goes here
%   Detailed explanation goes here

x = pos(1);
y = pos(2);

theta = ori;

a_e = ellipse(1);
b_e = ellipse(2);

beta = linspace(0,2*pi,res)';

r = a_e*b_e./(sqrt(b_e^2*cos(beta).*cos(beta)+a_e^2*sin(beta).*sin(beta)));

v = r .* cos(beta);
w = r .* sin(beta);

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

ic = zeros(2,length(beta));

for i=1:length(beta)
    
    ic(:,i) = [x;y] + R*[v(i);w(i)];
    
end

[Y_min, ind ] = min(ic(2,:));
beta_c = beta(ind);

end

