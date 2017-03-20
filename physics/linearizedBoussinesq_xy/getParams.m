function p = getParams(p)

    % Linearized Boussinesq equations
    p.nVars = 4;
    p.cn = p.f0/p.kappa;

    % Domain wave number
    k1 = 2*pi/p.Lx;
    l1 = 2*pi/p.Ly;

    % x-grid from 0 to L-dx, centered on L/2.
    dx		= p.Lx/p.nx; 
    p.xx	= 0:dx:p.Lx-dx; p.xx=p.xx-p.Lx/2;
    p.dx	= p.xx(2)-p.xx(1);

    dy		= p.Ly/p.ny; 
    p.yy	= (0:dy:p.Ly-dy)'; p.yy=p.yy-p.Ly/2;
    p.dy	= p.yy(2)-p.yy(1);

    % Fourier wavenumbers
    p.kk = k1*[0:p.nx/2 -p.nx/2+1:-1];
    p.ll = l1*[0:p.ny/2 -p.ny/2+1:-1]';

    % Mesh of k's and l's
    [p.KK, p.LL] = meshgrid(p.kk, p.ll);
    [p.XX, p.YY] = meshgrid(p.xx, p.yy);

    p.kay2 = p.KK.^2 + p.LL.^2;
    p.kay2(1, 1) = Inf;

    % Zero out high wavenumbers
    p.filt  = zeros(p.ny, p.nx);

    % Cut-off wavenumber.
    kkc = 0.65*max(abs(p.kk));
    llc = 0.65*max(abs(p.ll));
    p.filt( (p.KK/kkc).^2 + (p.LL/llc).^2 < 1 ) = 1;

    p.filt = p.filt(:, :, ones(1, p.nVars));

end
