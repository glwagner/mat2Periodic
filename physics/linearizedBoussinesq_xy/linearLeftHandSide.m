function mu = linearLeftHandSide(p, KK, LL)

% Initialize
mu = zeros(p.ny, p.nx, p.nVars);

% Exponent reflecting the order of the diffusivity operator
pow = floor(p.nuN/2);

% Diffusivity operator
D = ( p.nu^(1/pow)*KK.^2 + p.nu^(1/pow)*LL.^2 ).^pow;

% Hyperviscous linear operator for q
mu(:, :, 4) = D;
