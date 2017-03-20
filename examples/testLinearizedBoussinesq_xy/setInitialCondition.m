function sol = setInitialCondition(p)

    % Initial condition
    %[sol0, icType] = randomCosines(p);
    %[sol0, icType] = randomNumbers(p);

    rmsRo = 0.05;
    [sol0, icType] = randomEnergySpectrum(p, rmsRo);

    switch icType
        case 'Physical'
            sol = p.filt.*fft2(sol0);
        case 'Spectral'
            sol = p.filt.*sol0;
    end
end

function [sol0, icType] = randomCosines(p)

    icType = 'Physical';

    % Initialize solution
    sol0 = zeros(p.ny, p.nx, p.nVars);

    % Random cosine parameters
    kz = 2*pi/p.Lx*(2:4);
    px = 2*pi*[0, 1/2, 5/6];
    py = 2*pi*[1/3, 3/7, 6/7];
    aa = [0.7, 0.2, 0.1];

    sol0x = zeros(p.ny, p.nx);
    sol0y = zeros(p.ny, p.nx);
    for ii = 1:length(kz)
        sol0x = sol0x + aa(ii)*cos(kz(ii)*p.XX+px(ii));
        sol0y = sol0y + aa(ii)*cos(kz(ii)*p.YY+py(ii));
    end

    sol0(:, :, 4) = sol0x.*sol0y;
    sol0(:, :, 4) = sol0(:, :, 4) ./ max(max(sol0(:, :, 4)));

    % Horizontal uniform inertial oscillation
    sol0(:, :, 1) = ones(p.ny, p.nx); 

end


function [sol0, icType] = randomNumbers_uniformWave(p)

    icType = 'Physical';

    sol0 = zeros(p.ny, p.nx, p.nVars);

    % Random numbers
    sol0(:, :, 4) = rand(p.ny, p.nx); 

    % Horizontal uniform inertial oscillation
    sol0(:, :, 1) = ones(p.ny, p.nx); 

end

function [sol0, icType] = randomEnergySpectrum(p, rmsRo)

    icType = 'Physical';

    sol0 = zeros(p.ny, p.nx, p.nVars);

    % Non-dimensional isotropic wavenumber
    k  = sqrt(p.KK.^2+p.LL.^2)*p.Lx/(2*pi);

    % Energy spectral with peak at non-dim wavenumber 'kPeak'
    kPeak = 8;
    Ek = k.^4 ./ (1+k/kPeak).^16;

    % Initial vorticity
    qh = sqrt(Ek).*sqrt(p.KK.^2+p.LL.^2).*exp(1i*2*pi*rand(p.ny, p.nx));
    qh(1, 1) = 0;

    % Rescale q0 so Rossby number is Ro
    qrms = sqrt(mean(mean( real(ifft2(qh)).^2 )));
    sol0(:, :, 4) = p.f0*rmsRo*real(ifft2(qh))/qrms;

    % Horizontal uniform inertial oscillation
    sol0(:, :, 1) = ones(p.ny, p.nx); 

end
