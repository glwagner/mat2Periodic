function sol = setInitialCondition(p)

    % Initial condition
    [sol0, icType] = randomCosines(p);
    %[sol0, icType] = randomNumbers(p);
    %[sol0, icType] = randomEnergySpectrum(p)

    switch icType
        case 'Physical'
            sol = p.filt.*fft2(sol0);
        case 'Spectral'
            sol = p.filt.*sol0;
    end
end

function [sol0, icType] = randomCosines(p)

    icType = 'Physical';

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

    sol0 = sol0x.*sol0y;
    sol0 = sol0 ./ max(max(sol0));

end


function [sol0, icType] = randomNumbers(p)

    icType = 'Physical';

    % Random numbers
    sol0 = rand(p.ny, p.nx); 

end

function [sol0, icType] = randomEnergySpectrum(p)

    icType = 'Physical';

    % Non-dimensional isotropic wavenumber
    k  = sqrt(KK.^2+LL.^2)*p.Lx/(2*pi);

    % Energy spectral with peak at non-dim wavenumber 'kPeak'
    kPeak = 8;
    Ek = k.^4 ./ (1+k/kPeak).^16;

    % Initial vorticity
    qh = sqrt(Ek).*sqrt(KK.^2+LL.^2).*exp(1i*2*pi*rand(p.ny, p.nx));
    qh(1, 1) = 0;

    % Rescale q0 so rms eddy turnover time is "ep"
    ep = 0.1;
    qrms = sqrt(mean(mean( real(ifft2(qh)).^2 )));
    sol0 = ep*real(ifft2(qh))/qrms;

end


