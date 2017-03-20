function quickPlot(p, sol)

% Get vorticity
q = real(ifft2(sol));

% Define index vectors
nn = -p.nx/2+1:p.nx/2;
mm = -p.ny/2+1:p.ny/2;

figure(1), clf

ax(1) = subplot(1, 2, 1);
pcolor(p.xx, p.yy, q), shading flat

ax(1) = subplot(1, 2, 1);
pcolor(p.xx, p.yy, q), shading flat

ax(2) = subplot(1, 2, 2);
pcolor(nn, mm, fftshift(abs(sol))), shading flat

xlabel(ax(1), 'x'), ylabel(ax(1), 'y'), colorbar(ax(1), 'northoutside')
xlabel(ax(2), 'k'), ylabel(ax(2), 'l'), colorbar(ax(2), 'northoutside')

warning off
colormap(ax(1), cbrewer('div', 'RdBu', 64))
colormap(ax(2), flipud(cbrewer('seq', 'YlGnBu', 64)))

ax(2).XLim = round([-1 1]*1/3*p.nx);
ax(2).YLim = round([-1 1]*1/3*p.ny);

drawnow; pause(0.01)
