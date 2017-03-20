function mu = linearLeftHandSide(p, KK, LL)

% Exponent reflecting the order of the diffusivity operator
pow = floor(p.nuN/2);

% Diffusivity operator
D = ( p.nu^(1/pow)*KK.^2 + p.nu^(1/pow)*LL.^2 ).^pow;

% Viscous linear operator
mu = D;
