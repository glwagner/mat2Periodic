function diags = getDiagnostics(p, sol)

% Solution key
phih = sol(:, :, 1);
qh   = sol(:, :, 2);

phi0h = p.sol0(:, :, 1);
q0h   = p.sol0(:, :, 2);


diags.ft.name = 'Time in inertial periods';
diags.ft.units = [];
diags.ft.print = 1;

diags.ft.value = p.t * p.f/(2*pi);


diags.mKE.name = 'Mean kinetic energy';
diags.mKE.units = 'm^2/s^2';
diags.mKE.print = 0;

diags.mKE.value = 1/(p.nx*p.ny)*sum(sum( ...
                    1/2 * abs(qh).^2 ./ p.kay2 ...
                                       )); 


diags.dmKE.name = 'Change in mean kinetic energy';
diags.dmKE.units = [];
diags.dmKE.print = 1;

KE0 = 1/(p.nx*p.ny)*sum(sum( ...
                    1/2 * abs(q0h).^2 ./ p.kay2 ...
                                       )); 
diags.dmKE.value = diags.mKE.value / KE0;


diags.maxRo.name = 'Maximum Rossby number';
diags.maxRo.units = [];
diags.maxRo.print = 1;

diags.maxRo.value = max(max(abs(real(ifft2(qh))/p.f)));


diags.wKE.name = 'Wave kinetic energy';
diags.wKE.units = 'm^2/s^2';
diags.wKE.print = 0;

diags.wKE.value = 0.0;


diags.dwKE.name = 'Change in wave energy';
diags.dwKE.units = [];
diags.dwKE.print = 0;

diags.dwKE.value = 0.0;
