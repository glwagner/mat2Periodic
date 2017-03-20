function wavePlots(p,sol)

% limits for spectra plots
waveSpectraLimz = [-1 1]*24;
psiSpectraLimz  = [-1 1]*24;

% parameters
fs = 14;

% spectral matrices
[KK,LL]=meshgrid(p.kk,p.ll);
nn = -p.nx/2:p.nx/2-1;
mm = -p.ny/2:p.ny/2-1;

% vorticity in real space
q = real(ifft2(sol.qh));
qrms = sqrt(mean(mean( q.^2 )));
qmax = max(max( abs(q) ));

% kinetic energy spectra of q
psiSpec = abs(sol.qh);

% normalize
psiSpec = psiSpec / sum(sum(psiSpec));

% wave velocities in real space
uh = -1/(p.alpha*p.f0)*( ...
			1i*p.sig*bsxfun(@times,1i*p.kk,sol.ah) ...
			- p.f0*bsxfun(@times,1i*p.ll,sol.ah) );

vh = -1/(p.alpha*p.f0)*( ...
			1i*p.sig*bsxfun(@times,1i*p.ll,sol.ah) ...
			+ p.f0*bsxfun(@times,1i*p.kk,sol.ah) );

% transform to real space and add complex conjugate
u = ifft2(uh);
u = real(u + conj(u));

v = ifft2(vh);
v = real(v + conj(v));

% wave speed
sp = sqrt(u.^2+ v.^2);
sph = abs(uh).^2 + abs(vh).^2;

% pressure amplitude
a = ifft2(sol.ah);

% what to plot?
wavePlot = sp;
%wavePlot = real(a + conj(a));
waveSpec = abs(sol.ah);

% normalize
waveSpec = waveSpec / sum(sum(waveSpec));

% plot
hfig=figure(p.fig); clf, set(gcf,'DefaultTextInterpreter','latex')

ax(1)=subplot(2,2,1);
pcolor(p.xx*1e-3,p.yy*1e-3,q/p.f0), shading interp, axis square
%xlabel('$x$ (km)'), ylabel('$y$ (km)')
ylabel('$y$ (km)')

ax(2)=subplot(2,2,2); 
pcolor(p.xx*1e-3,p.yy*1e-3,wavePlot), shading interp, axis square
%xlabel('$x$ (km)'), ylabel('$y$ (km)')
ylabel('$y$ (km)')

ax(3)=subplot(2,2,3); 
pcolor(nn,mm,fftshift(psiSpec)), shading flat, axis square
xlabel('$k$'), ylabel('$\ell$')
ax(3).XLim = waveSpectraLimz;
ax(3).YLim = waveSpectraLimz;

ax(4)=subplot(2,2,4); 
pcolor(nn,mm,fftshift(waveSpec)), shading flat, axis square
xlabel('$k$'), ylabel('$\ell$')
ax(4).XLim = psiSpectraLimz;
ax(4).YLim = psiSpectraLimz;

% axis
ax(1).XTickLabel = [];
ax(2).XTickLabel = [];

% make colorbars
hc(1)=colorbar(ax(1),'northoutside');
hc(2)=colorbar(ax(2),'northoutside');

% color axis
caxis(ax(1),[-1 1]*qrms/p.f0*2)
caxis(ax(2),[0.01 2])
caxis(ax(3),caxis(ax(3))/2)
caxis(ax(4),caxis(ax(4))/4)

% change colormap
warning off
colormap(ax(1),polarmap(jet));
colormap(ax(2),flipud(cbrewer('seq', 'YlGnBu', 64)));
colormap(ax(3),flipud(cbrewer('seq', 'YlGnBu', 64)));
colormap(ax(4),flipud(cbrewer('seq', 'YlGnBu', 64)));

% text labels
dxtxt = 0.01;
dytxt = 0.04;
txt3 = ['$t = ' num2str(round(p.sig*p.t/(2*pi)),'%04d') '$'];
htxt1 = text(dxtxt,dytxt,'$\triangle \psi /f_0$', ...
				'Units','Normalized','Parent',ax(1), ...
				'FontSize',fs);
%htxt2 = text(dxtxt,dytxt,'$A+A^*$', ...
htxt2 = text(dxtxt,dytxt,'$\sqrt{\tilde{u}^2+\tilde{v}^2}$', ...
				'Units','Normalized','Parent',ax(2), ...
				'FontSize',fs,'Color','w');
htxt3 = text(1-dxtxt,dytxt,txt3, ...
				'Units','Normalized','Parent',ax(2), ...
				'FontSize',fs,'HorizontalAlignment','right','Color','w');
%htxt4 = text(dxtxt,dytxt,'$\left ( k^2 + \ell^2 \right ) | \hat \psi |^2$', ...
htxt4 = text(dxtxt,dytxt,'$| \hat q |$', ...
				'Units','Normalized','Parent',ax(3), ...
				'FontSize',fs,'Color','w');
htxt5 = text(dxtxt,dytxt,'$|\hat A|$', ...
				'Units','Normalized','Parent',ax(4), ...
				'FontSize',fs,'Color','w');
% pretty
ax(1).FontSize=fs;
ax(2).FontSize=fs;
ax(3).FontSize=fs;
ax(4).FontSize=fs;

ax(2).YAxisLocation = 'right';
ax(4).YAxisLocation = 'right';

% readjust for colorbars
ax(1).Position(4) = ax(3).Position(4);
ax(2).Position(4) = ax(3).Position(4);

% squeeze horizontally
xsqueeze = 0.1;
ax(2).Position(1) = ax(2).Position(1)-xsqueeze;
ax(4).Position(1) = ax(4).Position(1)-xsqueeze;

% shift vertically
vshift = 0.1;
ax(1).Position(2) = ax(1).Position(2)-vshift;
ax(2).Position(2) = ax(2).Position(2)-vshift;

%{
% shift colorbars slightly
shift = 0.01;
hc(1).Position(1) = hc(1).Position(1)+shift;
hc(2).Position(1) = hc(2).Position(1)+shift;
%}

% tick label interpreters
ax(1).TickLabelInterpreter = 'latex';
ax(2).TickLabelInterpreter = 'latex';
ax(3).TickLabelInterpreter = 'latex';
ax(4).TickLabelInterpreter = 'latex';
hc(1).TickLabelInterpreter = 'latex';
hc(2).TickLabelInterpreter = 'latex';

% more analysis

% rms vorticity
q = real(ifft2(sol.qh));
qrms = sqrt(mean(mean(q.^2)));
qmax = max(max(abs(q)));

% energy containing length scale
psih = -sol.qh./p.kay; psih(1,1) = 0;
numer = trapz(p.kk,trapz(p.ll,p.kay.*abs(psih).^2));
denom = trapz(p.kk,trapz(p.ll,p.kay.^(3/2).*abs(psih).^2));
Lbar = numer/denom;

disp(sprintf('\nRo_rms = %5.3f, Ro_max = %5.3f, Lbar = %4.1f km \n', ...
				qrms/p.f0,qmax/p.f0,Lbar*1e-3))
