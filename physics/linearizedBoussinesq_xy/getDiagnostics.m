function diags = getDiagnostics(p, sol)

    % Solution key
    uh = sol(:, :, 1);
    vh = sol(:, :, 2);
    hh = sol(:, :, 3);
    qh = sol(:, :, 4);

    u0h = p.sol0(:, :, 1);
    v0h = p.sol0(:, :, 2);
    h0h = p.sol0(:, :, 3);
    q0h = p.sol0(:, :, 4);


    diags.ft.name = 'Time in inertial periods';
    diags.ft.units = [];
    diags.ft.print = 1;

    diags.ft.value = p.t * p.f0/(2*pi);


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

    diags.wKE.name = 'Wave kinetic energy';
    diags.wKE.units = 'm^2/s^2';
    diags.wKE.print = 0;

    diags.wKE.value = 1/(p.nx*p.ny)*sum(sum( ...
                        1/2 * (abs(uh).^2 + abs(vh).^2) ...
                                           ));

    diags.wPE.name = 'Wave potential energy';
    diags.wPE.units = 'm^2/s^2';
    diags.wPE.print = 0;

    diags.wPE.value = 1/(p.nx*p.ny)*sum(sum( ...
                        1/2 * abs(hh).^2 / p.cn^2 ...
                                           ));

    diags.wE.name = 'Total wave energy';
    diags.wE.units = 'm^2/s^2';
    diags.wE.print = 0;

    diags.wE.value = diags.wKE.value + diags.wPE.value;


    diags.wKE0.name = 'Initial wave kinetic energy';
    diags.wKE0.units = 'm^2/s^2';
    diags.wKE0.print = 0;

    diags.wKE0.value = 1/(p.nx*p.ny)*sum(sum( ...
                        1/2 * (abs(u0h).^2 + abs(v0h).^2) ...
                                           ));


    diags.wPE0.name = 'Initial wave potential energy';
    diags.wPE0.units = 'm^2/s^2';
    diags.wPE0.print = 0;

    diags.wPE0.value = 1/(p.nx*p.ny)*sum(sum( ...
                        1/2 * abs(h0h).^2 / p.cn^2 ...
                                           ));


    diags.dwE.name = 'Change in total wave energy';
    diags.dwE.units = [];
    diags.dwE.print = 1;

    diags.dwE.value = diags.wE.value/(diags.wKE0.value + diags.wPE0.value);

end
