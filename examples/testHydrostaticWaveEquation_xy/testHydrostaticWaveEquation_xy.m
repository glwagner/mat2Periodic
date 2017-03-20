clear all

p.home = '/Users/glwagner/Numerics/mat2Periodic';
p.physics = 'hydrostaticWaveEquation_xy';
p.name = 'testHydrostaticWaveEquation_xy';

p.nx = 128;
p.Lx = 1e6;
p.nSteps = 1e3;

% Physical parameters: inertial freq, kappa, viscosity
p.f0 = 1e-4;
p.alpha = 3;
p.sigma = sqrt(1+p.alpha)*p.f0;
% Recall that k = sqrt(alpha)*kappa
p.kappa = (32*pi/p.Lx) / sqrt(p.alpha);

p.nuQ = 1e6;
p.nuA = 1e2;
p.nuNQ = 4;
p.nuNA = 4;

p.dt = 5e-2 * (2*pi/p.sigma);

% Step-interval between display.
p.dnPrint = 2e2;
p.dnPlot = p.dnPrint;
p.dnSave = p.dnPrint;

%p.timeStepper = 'RK4';
p.timeStepper = 'ETDRK4';

% Square domain
p.ny = p.nx;
p.Ly = p.Lx;

% Manage code in the sadly-MATLABian-necessary way.
codeDir = sprintf('%s/code', pwd);
sourceDir = sprintf('%s/source', p.home);
physicsDir = sprintf('%s/physics/%s', p.home, p.physics);

if ~exist(codeDir), mkdir(codeDir)
else eval(sprintf('!rm %s/*', codeDir))
end

addpath(codeDir)

% Copy code and mat-files into code directory
eval(sprintf('!cp %s/*.m %s/', sourceDir, codeDir)) 
eval(sprintf('!cp %s/*.m %s/', physicsDir, codeDir)) 
eval(sprintf('!cp %s/*.m* %s/', pwd, codeDir)) 

% Initialize and run the model
[p, sol] = initializeModel(p);
[p, sol] = runSpectralModel(p, sol, p.nSteps);

% Plot solution at the end.
quickPlot(p, sol);
