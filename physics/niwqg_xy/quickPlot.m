function quickPlot(p, sol)

% Solution key
Ah = sol(:, :, 1);
qh = sol(:, :, 2);

fontSize = 16;

q = real(ifft2(qh));
A = ifft2(Ah);

figure(1), clf

ax(1) = subplot(121);
pcolor(p.xx, p.yy, q/p.f), shading flat

ax(2) = subplot(122);
pcolor(p.xx, p.yy, abs(A)), shading flat

colormap(ax(1), flipud(cbrewer('div', 'RdBu', 64)));
colormap(ax(2), flipud(cbrewer('seq', 'YlGnBu', 64)));

hc(1) = colorbar(ax(1), 'northoutside');
hc(2) = colorbar(ax(2), 'northoutside');

pause(0.1)
