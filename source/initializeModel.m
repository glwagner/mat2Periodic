function [p, sol] = initializeModel(p)

    fprintf('Initializing... '), t1=tic;

    if isfield(p, 'mat2PeriodicDir')

        % Manage code in the sadly-MATLABian-necessary way.
        p.codeDir = sprintf('%s/code', p.homeDir);
        p.physicsDir = sprintf('%s/physics/%s', p.mat2PeriodicDir, p.physics);

        if ~exist(p.codeDir), mkdir(p.codeDir), end
        addpath(p.codeDir)

        % Copy code and mat-files into code directory
        eval(sprintf('!cp %s/*.m %s/', p.sourceDir, p.codeDir)) 
        eval(sprintf('!cp %s/*.m %s/', p.physicsDir, p.codeDir)) 
        eval(sprintf('!cp %s/*.m* %s/', pwd, p.codeDir)) 

    end

    % Get problem parameters.
    p = getParams(p);

    % Set initial condition.
    sol = setInitialCondition(p); 
    p.sol0 = sol;

    % Initialze linear left hand side.
    p.mu = linearLeftHandSide(p, p.KK, p.LL);

    % Initialize the time-stepper.
    p = initializeTimeStepper(p, p.mu);

    fprintf('\\m/ (t = %6.3f s).\n', toc(t1))

end
