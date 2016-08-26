
% u is the contact-space velocity ([tan1 tan2 normal])
% mu is the coefficient of friction
% epsilon is the coefficient of restitution
function [v_plus, z] = smith(M, n, s, v, ha, mu, epsilon)

  % get the normal velocity
  nv = n'*(v + ha);

  % if nv is non-negative, quit
  if (nv >= 0)
    z = [0 0 s'*(v + ha)];
    v_plus = v;
    return;
  end

  % get the normal and tangent velocity
  nsv = [n s]'*(v + ha);

  % compute the contact space inertia matrix
  K = [n s]'*(M \ [n s]);

  % solve the nonlinear systme of equations for pn and vft
  x0 = [0; 0];
  x = fsolve(@(x) smith_fun(x, K, nsv, mu, epsilon), x0);

  % compute the collisional impulse
  pn = x(1);
  vft = x(2);
  vit = nsv(2);
  j = pn*[1; -mu*(vit*abs(vit)+vft*abs(vft))/(vit^2 + vft^2)];

  % compute the new velocity
  z = j;
  v_plus = M \ [n s]*z + (v + ha);

  % finish setting up z
  z(3) = 0;
  if (z(2) < 0)
    z(3) = -z(2);
    z(2) = 0;
  end

  % compute the new velocity
  z(4) = abs(s'* v_plus);

