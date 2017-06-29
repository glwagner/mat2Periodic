function sol = setInitialCondition(p)

    % Initialize solution
    sol0 = zeros(p.ny, p.nx, p.nVars);

    % Initial condition
    %sol0 = gaussianVortex(p, sol0, p.Lx/20, 0.05);
    sol0 = vortexCondensate(p, sol0, 0.1);
    %sol0 = randomEnergySpectrum(p, sol0, 32, 0.02);
    sol0 = uniformInertialOscillation(p, sol0);
    %sol0 = planeWave(1, p, sol0);
    
    % Transform into spectral space.
    sol = p.filt.*fft2(sol0);

end

function sol0 = uniformInertialOscillation(p, sol0)
    sol0(:, :, 1) = ones(p.ny, p.nx);
end

function sol0 = planeWave(n, p, sol0)
    k = 2*pi*n/p.Lx;
    sol0(:, :, 1) = exp(1i*k*p.XX);
end

function sol0 = gaussianVortex(p, sol0, R, Ro)
    sol0(:, :, 2) = -p.f*Ro*exp(-(p.XX.^2+p.YY.^2)/(2*R^2));
end

function sol0 = vortexCondensate(p, sol0, Ro)
    k = 2*pi/p.Lx;
    sol0(:, :, 2) = p.f*Ro*(sin(k*p.XX) + sin(k*p.YY));
end

function sol0 = randomNumbers(p, sol0, rmsRo)
    qh = rand(p.ny, p.nx); 
    q = real(ifft2(qh));
    sol(:, :, 2) = p.f * rmsRo * qh / sqrt(mean(mean(q.^2)));
end

function sol0 = randomEnergySpectrum(p, sol0, kPeak, rmsRo)

    % Non-dimensional isotropic wavenumber
    k  = sqrt(p.KK.^2+p.LL.^2)*p.Lx/(2*pi);

    % Energy spectral with peak at non-dim wavenumber 'kPeak'
    Ek = k.^4 ./ (1+k/kPeak).^16;

    % Initial vorticity
    qh = sqrt(Ek).*sqrt(p.KK.^2+p.LL.^2).*exp(1i*2*pi*rand(p.ny, p.nx));
    qh(1, 1) = 0;

    % Rescale q0 so Rossby number is Ro
    qrms = sqrt(mean(mean( real(ifft2(qh)).^2 )));
    sol0(:, :, 2) = p.f*rmsRo*real(ifft2(qh))/qrms;

end
