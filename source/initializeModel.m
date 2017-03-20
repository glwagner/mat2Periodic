function [p, sol] = initializeModel(p)

    fprintf('Initializing... '), t1=tic;

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
