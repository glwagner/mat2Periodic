function diags = getDiagnostics(p, sol)

    % Solution key
    Ah = sol(:, :, 1);
    qh = sol(:, :, 2);

    A0h = p.sol0(:, :, 1);
    q0h = p.sol0(:, :, 2);


    diags.sigt.name = 'Time in wave periods';
    diags.sigt.units = [];
    diags.sigt.print = 1;

    diags.sigt.value = p.t * p.sigma/(2*pi);


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


    diags.A.name = 'Wave action';
    diags.A.units = 'm^2/s^2';
    diags.A.print = 0;

    diags.A.value = 1/(p.nx*p.ny) * sum(sum( ...
                        1/(2*p.alpha*p.sigma) * ( ...
                            abs(p.KK.*Ah).^2 + abs(p.LL.*Ah).^2 ...
                          + (4+3*p.alpha)*p.kappa^2*abs(Ah).^2 ) ...
                                           ));


    diags.wE.name = 'Wave energy';
    diags.wE.units = 'm^2/s^2';
    diags.wE.print = 0;

    diags.wE.value = 1/(p.nx*p.ny) * sum(sum( ...
                        (p.alpha+2)/p.alpha^2 ...
                         * (abs(p.KK.*Ah).^2 + abs(p.LL.*Ah).^2) ...
                          + p.kappa^2*abs(Ah).^2 ...
                                            ));


    diags.dwE.name = 'Change in wave energy';
    diags.dwE.units = [];
    diags.dwE.print = 1;

    wE0 = 1/(p.nx*p.ny) * sum(sum( ...
                        (p.alpha+2)/p.alpha^2 ...
                         * (abs(p.KK.*A0h).^2 + abs(p.LL.*A0h).^2) ...
                          + p.kappa^2*abs(A0h).^2 ...
                                            ));

    diags.dwE.value = diags.wE.value / wE0;




    diags.dA.name = 'Change in wave action';
    diags.dA.units = [];
    diags.dA.print = 1;

    A0 = 1/(p.nx*p.ny) * sum(sum( ...
                        1/(2*p.alpha*p.sigma) * ( ...
                            abs(p.KK.*A0h).^2 + abs(p.LL.*A0h).^2 ...
                          + (4+3*p.alpha)*p.kappa^2*abs(A0h).^2 ) ...
                                ));

    diags.dA.value = diags.A.value / A0;

end
