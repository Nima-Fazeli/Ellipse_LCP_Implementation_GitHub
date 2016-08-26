function residual = smith_fun(x, K, u, mu, epsilon)

  % get pn and vft
  pn = x(1);
  vft = x(2);

  % get vin and vit
  vin = u(1);
  vit = u(2);

  % evaluate
  residual = K * [-(1+epsilon)*vin; vft - vit] - pn * [1; -mu*(vit*abs(vit)+vft*abs(vft))/(vit^2 + vft^2)];

