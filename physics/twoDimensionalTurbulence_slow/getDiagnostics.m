function diags = getDiagnostics(p, sol)

    % Time
    diags.t.name = 'time';
    diags.t.units = 's';
    diags.t.print = 1;
    
    diags.t.value = p.t;

    % Kinetic energy
    diags.KE.name = 'kinetic energy';
    diags.KE.units = 'm^2/s^2';
    diags.KE.print = 0;

    psih = sol./p.kay2;
    diags.KE.value = 1/(p.nx*p.ny) * 1/2 * sum(sum( ...
                      (p.KK.^2+p.LL.^2).*abs(psih).^2));

    % Kinetic energy
    diags.dKE.name = 'normalized kinetic energy';
    diags.dKE.units = [];
    diags.dKE.print = 1;

    psi0h = p.sol0 ./ p.kay2;
    diags.dKE.value = diags.KE.value / ...
                      (1/(p.nx*p.ny) * 1/2 * sum(sum( ...
                        (p.KK.^2+p.LL.^2).*abs(psi0h).^2)));


end
