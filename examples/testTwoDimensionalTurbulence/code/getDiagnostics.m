function diags = getDiagnostics(p, sol)

    % Kinetic energy
    diags.KE.name = 'kinetic energy';
    diags.KE.units = 'm^2/s^2';
    diags.KE.print = 1;

    psih = sol./p.kay2;
    diags.KE.value = 1/(p.nx*p.ny) * 1/2 * sum(sum( ...
                      (p.KK.^2+p.LL.^2).*abs(psih).^2));

end
