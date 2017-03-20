function p = initializeTimeStepper(p, mu)

    switch p.timeStepper
        case 'ETDRK4'
            p.ETDRK4Coeffs = getETDRK4Coeffs(p, mu);
    end

end

function coeffs = getETDRK4Coeffs(p, mu) 

    % Take mean over a circle with ncirc points
    ncirc	= 16;
    vcirc	= shiftdim(1:ncirc, -2);

    % Unit circle from e^{pi*i/n} to e^{2*pi*i*(1-1/(2n))}
    circ 	= exp(2*pi*1i*( vcirc-1/2 )/ncirc);

    % Draw unit circle around each mu_a and mu_b
    zc 	= -p.dt*mu(:, :, :, ones(ncirc,1)) ...
                + circ(ones(p.ny,1), ones(1,p.nx), ones(1,p.nVars), :);

    % Prefactor for substeps.  note that L = -mu
    zeta 	= p.dt*mean(( exp(zc/2)-1 )./zc , 4);

    % Prefactors for final step.
    alph 	= p.dt*mean(( -4 - zc + exp(zc).*(4-3*zc+zc.^2)	  )./zc.^3 , 4);
    beta	= p.dt*mean((  2 + zc + exp(zc).*(-2+zc)          )./zc.^3 , 4);
    gamm	= p.dt*mean(( -4 - 3*zc - zc.^2 + exp(zc).*(4-zc) )./zc.^3 , 4);

    % Store 
    coeffs.alph = alph;
    coeffs.beta = beta;
    coeffs.gamm = gamm;
    coeffs.zeta = zeta;
    coeffs.expMuDt = exp(-p.dt*mu/2);

end
