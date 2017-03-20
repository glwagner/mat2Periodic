clear all

p.physics = 'twoDimensionalTurbulence';
p.dt = 1e-1;
p.nx = 384;     
p.Lx = 2*pi;    
p.nSteps = 1e3;

% Viscosity and viscous order
p.nu = 1e-4;
p.nuN = 2;

% Step-interval between display.
p.dnPrint = 1e2;
p.dnPlot = p.dnPrint;
p.dnSave = 0;

p.timeStepper = 'RK4';
%p.timeStepper = 'ETDRK4';

% Square domain
p.ny = p.nx;
p.Ly = p.Lx;

% Manage code in the sadly-MATLABian-necessary way.
codeDir = sprintf('%s/code', pwd);
sourceDir = '../../source';
physicsDir = sprintf('../../physics/%s', p.physics);

if ~exist(codeDir), mkdir(codeDir), end
addpath(codeDir)

% Copy code and mat-files into code directory
eval(sprintf('!cp %s/*.m %s/', sourceDir, codeDir)) 
eval(sprintf('!cp %s/*.m %s/', physicsDir, codeDir)) 
eval(sprintf('!cp %s/*.m* %s/', pwd, codeDir)) 

% Initialize the model
[p, sol] = initializeModel(p);

% Run the model
[p, sol] = runSpectralModel(p, sol, p.nSteps);

% Plot solution at the end.
quickPlot(p, sol);
