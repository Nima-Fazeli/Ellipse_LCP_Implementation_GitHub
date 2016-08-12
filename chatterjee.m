
% u is the contact-space velocity ([tan1 tan2 normal])
% mu is the coefficient of friction
% epsilon is the coefficient of restitution
function [v_plus, z] = chatterjee(M, n, s, v, ha, mu, epsilon, epsilon_t)

  % get the normal velocity
  nv = n'*(v + ha);

  % if nv is non-negative, quit
  if (nv >= 0)
    z = [0 0 s'*(v + ha)];
    v_plus = v;
    return;
  end

  % compute the contact space inertia matrix along the contact normal
  nmnT = n'*(M \ n);

  % compute the perfectly plastic frictionless impulse
  z_I = nmnT \ -nv;
  z_I(2,1) = 0;

  % get the normal and tangent velocity
  nsv = [n s]'*(v + ha);

  % compute the contact space inertia matrix
  K = [n s]'*(M \ [n s]);

  % compute the perfectly plastic, full friction impulse
  z_II = K \ -nsv;

  % compute the perfectly plastic frictionless impulse
  z_hat = (1+epsilon)*z_I + (1+epsilon_t)*(z_II - z_I);

  % friction cone check
  if (abs(z_hat(2)) > mu*abs(z_hat(1)))
    kappa = (-z_I(2) - epsilon*z_I(2) + mu*z_I(1) + epsilon*mu*z_I(1))/(-mu*z_II(1) + z_II(2) + mu*z_I(1) - z_I(2)); 
  else
    kappa = 1 + epsilon_t;
  end

  % set the collisional impulse
  z = (1+epsilon)*z_I + kappa*(z_II - z_I);
  mu*z(1) - abs(z(2))

  % compute the new velocity
  v_plus = M \ [n s]*z + (v + ha);
  z(4) = abs(s'* v_plus);

