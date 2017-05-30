clear all

p.physics       = 'hydrostaticWaveEquation_xy';
p.Lx            = 1e6;
p.nSteps        = 2e2;
p.f0            = 1e-4;
p.alpha         = 3;
p.sigma         = sqrt(1+p.alpha)*p.f0;
p.kappa         = (32*pi/p.Lx) / sqrt(p.alpha);
p.nuQ           = 1e6;
p.nuA           = 1e2;
p.nuNQ          = 4;
p.nuNA          = 4;
p.dt            = 5e-2 * (2*pi/p.sigma);
p.dnPrint       = 0;
p.dnPlot        = p.dnPrint;
p.dnSave        = p.dnPrint;
p.timeStepper   = 'ETDRK4';

p.homeDir = pwd;
p.mat2PeriodicDir = [pwd '/../..'];
p.sourceDir = [pwd '/../../source'];
addpath(p.sourceDir)

times = []; nxs = [];
for nx = 2.^[7, 8, 9, 10]

    p.nx = nx;     
    p.ny = p.nx;
    p.Ly = p.Lx;

    % Initialize the model
    [p, sol] = initializeModel(p);

    % Run the model
    t1 = tic;
    [p, sol] = runSpectralModel(p, sol, p.nSteps);

    times = [times toc(t1)];
    nxs = [nxs nx];

end

for i = 1:length(times)
    fprintf('\nnx: %d, time: %.6f secs', nxs(i), times(i));
end
fprintf('\n')
