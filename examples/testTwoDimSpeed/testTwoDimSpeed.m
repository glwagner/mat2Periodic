clear all

p.physics       = 'twoDimensionalTurbulence';
p.dt            = 1e-1;
p.Lx            = 2*pi;    
p.nu            = 1e-4;
p.nuN           = 2;
p.dnPrint       = 0;
p.dnPlot        = p.dnPrint;
p.dnSave        = 0;
p.timeStepper   = 'RK4';
p.nSteps        = 2e2;

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
