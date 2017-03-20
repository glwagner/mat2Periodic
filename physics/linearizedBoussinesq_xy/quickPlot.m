function quickPlot(p, sol)

    % Solution key
    uh = sol(:, :, 1);
    vh = sol(:, :, 2);
    hh = sol(:, :, 3);
    qh = sol(:, :, 4);

    fontSize = 16;

    % Get fields
    u = real(ifft2(uh));
    v = real(ifft2(vh));
    q = real(ifft2(qh));
    sp = 1/2*sqrt(u.^2 + v.^2);

    figure(1), clf
    set(gcf, 'DefaultTextInterpreter', 'latex')

    ax(1) = subplot(1, 2, 1);
    pcolor(1e-3*p.xx, 1e-3*p.yy, q), shading flat
    axis square

    ax(2) = subplot(1, 2, 2);
    pcolor(1e-3*p.xx, 1e-3*p.yy, sp/max(max(sp))), shading flat
    axis square

    xlabel(ax(1), '$x$')
    ylabel(ax(1), '$y$')
    cb(1) = colorbar(ax(1), 'northoutside');

    xlabel(ax(2), '$x$')
    ylabel(ax(2), '$y$')
    cb(2) = colorbar(ax(2), 'northoutside');

    ax(1).FontSize = fontSize;
    ax(2).FontSize = fontSize;

    ax(2).YAxisLocation = 'right';

    ax(1).TickLabelInterpreter = 'latex';
    ax(2).TickLabelInterpreter = 'latex';
    cb(1).TickLabelInterpreter = 'latex';
    cb(2).TickLabelInterpreter = 'latex';

    warning off
    colormap(ax(1), flipud(cbrewer('div', 'RdBu', 64)))
    colormap(ax(2), flipud(cbrewer('seq', 'YlGnBu', 64)))

    drawnow; pause(0.1)

    % Squeeze
    xSqueeze = 0.02;
    ax(1).Position(1) = ax(1).Position(1) + xSqueeze;
    ax(2).Position(1) = ax(2).Position(1) - xSqueeze;

    drawnow; pause(0.1)

end
