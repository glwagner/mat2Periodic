function mu = linearLeftHandSide(p, KK, LL)

% Initialize
mu = zeros(p.ny, p.nx, p.nVars);

% Exponent reflecting the order of the diffusivity operator
powQ = p.nuNQ/2;
powA = p.nuNA/2;

% Diffusivity operator
DA = ( p.nuA^(1/powA)*KK.^2 + p.nuA^(1/powA)*LL.^2 ).^powA;
DQ = ( p.nuQ^(1/powQ)*KK.^2 + p.nuQ^(1/powQ)*LL.^2 ).^powQ;

% Dispersion
Disp = 1i*p.f*(p.KK.^2+p.LL.^2)./(2*p.kappa^2) .* p.invL;

% Linear dispersion for A and hyperviscosity for A and q
mu(:, :, 1) = DA + Disp;
mu(:, :, 2) = DQ;
