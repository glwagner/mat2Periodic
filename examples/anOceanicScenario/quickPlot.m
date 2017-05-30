function quickPlot(p, sol)

    % Solution key
    Ah = sol(:, :, 1);
    qh = sol(:, :, 2);

    fontSize = 16;

    nn = -p.nx/2:p.nx/2-1;
    mm = -p.ny/2:p.ny/2-1;

    % Get fields
    q = real(ifft2(qh));
    [u, v] = getVelocities(p, Ah);
    sp = sqrt(u.^2 + v.^2);

    figure(1), clf
    set(gcf, 'DefaultTextInterpreter', 'latex')

    ax(1) = subplot(2, 2, 1);
    pcolor(1e-3*p.xx, 1e-3*p.yy, q/p.f0), shading flat
    axis square

    ax(2) = subplot(2, 2, 2);
    pcolor(1e-3*p.xx, 1e-3*p.yy, sp), shading flat
    axis square

    ax(3) = subplot(2, 2, 3);
    pcolor(nn, mm, fftshift(abs(qh))), shading flat
    axis square

    ax(4) = subplot(2, 2, 4);
    pcolor(nn, mm, fftshift(abs(Ah))), shading flat
    axis square

    drawnow; pause(0.1)

    ax(3).XLim = [-1 1]*32;
    ax(3).YLim = ax(3).XLim;

    ax(4).XLim = ax(3).XLim;
    ax(4).YLim = ax(3).XLim;

    cb(1) = colorbar(ax(1), 'northoutside');
    cb(2) = colorbar(ax(2), 'northoutside');

    ax(1).XTickLabels = [];
    ax(2).XTickLabels = [];

    ax(1).CLim = [-1 1]*0.1;
    ax(2).CLim = [0 2];
    ax(3).CLim = ax(3).CLim/2;
    ax(4).CLim = ax(4).CLim/2;

    xlabel(ax(3), '$k$')
    xlabel(ax(4), '$k$')

    ylabel(ax(1), '$y$')
    ylabel(ax(2), '$y$')
    ylabel(ax(3), '$\ell$')
    ylabel(ax(4), '$\ell$')

    ax(1).FontSize = fontSize;
    ax(2).FontSize = fontSize;
    ax(3).FontSize = fontSize;
    ax(4).FontSize = fontSize;

    ax(2).YAxisLocation = 'right';
    ax(4).YAxisLocation = 'right';

    ax(1).TickLabelInterpreter = 'latex';
    ax(2).TickLabelInterpreter = 'latex';
    ax(3).TickLabelInterpreter = 'latex';
    ax(4).TickLabelInterpreter = 'latex';

    cb(1).TickLabelInterpreter = 'latex';
    cb(2).TickLabelInterpreter = 'latex';

    warning off
    colormap(ax(1), flipud(cbrewer('div', 'RdBu', 64)))
    colormap(ax(2), flipud(cbrewer('seq', 'YlGnBu', 64)))
    colormap(ax(3), flipud(cbrewer('seq', 'YlGnBu', 64)))
    colormap(ax(4), flipud(cbrewer('seq', 'YlGnBu', 64)))

    drawnow; pause(0.1)

    % Restore dimensions
    ax(1).Position(1) = ax(3).Position(1);
    ax(2).Position(1) = ax(4).Position(1);

    ax(1).Position(3) = ax(3).Position(3);
    ax(1).Position(4) = ax(3).Position(4);

    ax(2).Position(3) = ax(3).Position(3);
    ax(2).Position(4) = ax(3).Position(4);

    % Squeeze
    xSqueeze = 0.06;
    ax(1).Position(1) = ax(1).Position(1) + xSqueeze;
    ax(3).Position(1) = ax(3).Position(1) + xSqueeze;
    ax(2).Position(1) = ax(2).Position(1) - xSqueeze;
    ax(4).Position(1) = ax(4).Position(1) - xSqueeze;

    ySqueeze = 0.04;
    ax(1).Position(2) = ax(1).Position(2) - ySqueeze;
    ax(2).Position(2) = ax(2).Position(2) - ySqueeze;
    ax(3).Position(2) = ax(3).Position(2) + ySqueeze;
    ax(4).Position(2) = ax(4).Position(2) + ySqueeze;

    drawnow; pause(0.1)

end
