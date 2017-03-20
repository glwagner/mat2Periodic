function mu = linearLeftHandSide(p, KK, LL)

% Initialize
mu = zeros(p.ny, p.nx, p.nVars);

% Exponent reflecting the order of the diffusivity operator
powQ = p.nuNQ/2;
powA = p.nuNA/2;

% Diffusivity operator
DQ = ( p.nuQ^(1/powQ)*KK.^2 + p.nuQ^(1/powQ)*LL.^2 ).^powQ;
DA = ( p.nuA^(1/powA)*KK.^2 + p.nuA^(1/powA)*LL.^2 ).^powA;

% Dispersion
Disp = 1i*p.alpha*p.sigma*p.invE.*(p.alpha*p.kappa^2 - p.KK.^2 - p.LL.^2);

% Linear dispersion for A and hyperviscosity for A and q
mu(:, :, 1) = DA + Disp;
mu(:, :, 2) = DQ;
