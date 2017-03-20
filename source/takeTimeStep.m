function sol = takeTimeStep(p, mu, sol)

    % Take a time-step using specified time-stepper. 
    % Note that linearly implicit time-steppers 
    % do not require mu, while linearly explicit steppers do.

    % Time-steppers:
    %   ETDRK4  : 4th-order Runge-kutta exponential time-marching. 
    %   RK4     : Fully explicit 4th-order Runge-Kutta time-marching.
    %   RKW3    : Low-storage 3rd-order Runge-Kutta time-marching.

    switch p.timeStepper
        case 'ETDRK4'
            sol = etdrk4TimeStep(p, sol);
        case 'RK4'
            sol = rk4TimeStep(p, mu, sol);
        case 'RKW3'
            sol = rkw3TimeStep(p, mu, sol);
    end

end

function sol = rkw3TimeStep(p, mu, sol)

    b1 = p.RKW3Coeffs.b1;
    b2 = p.RKW3Coeffs.b2;
    b3 = p.RKW3Coeffs.b3;

    c2 = p.RKW3Coeffs.c2;
    c3 = p.RKW3Coeffs.c3;

    a21 = p.RKW3Coeffs.a21;
    a32 = p.RKW3Coeffs.a32;

    RHS0 = rightHandSide(p, sol) - mu.*sol;

    p.t = p.t + c2*p.dt;
    intSol = a21*p.dt*RHS0;
    RHS1 = rightHandSide(p, intSol) - mu.*intSol;

    % Reverse previous step in p.t and then compute substep.
    p.t = p.t - c2*p.dt + c3*p.dt;
    intSol = b1*p.dt*RHS0 + a32*p.dt*RHS1;
    RHS2 = rightHandSide(p, intSol);

    % Reverse previous step in p.t and step forward.
    sol = sol + p.dt*( b1*RHS0 + b2*RHS1 + b3*RHS2 );

end
    

function sol = rk4TimeStep(p, mu, sol)

    % Take a 4th-order Runge-Kutta time-step.
    % Recall that mu.*sol occupies the left side of the 
    % differential equation.

    RHS0 = rightHandSide(p, sol) - mu.*sol;

    p.t = p.t + p.dt/2;
    intSol = sol + p.dt/2*RHS0;
    RHS1 = rightHandSide(p, intSol) - mu.*intSol;

    intSol = sol + p.dt/2*RHS1;
    RHS2 = rightHandSide(p, intSol) - mu.*intSol;

    p.t = p.t + p.dt/2;
    intSol = sol + p.dt*RHS2;
    RHS3 = rightHandSide(p, intSol) - mu.*intSol;

    sol = sol + p.dt*( 1/6*RHS0 + 1/3*RHS1 + 1/3*RHS2 + 1/6*RHS3 );

end

function sol = etdrk4TimeStep(p, sol);

    % Convenient views.
    alph = p.ETDRK4Coeffs.alph;
    beta = p.ETDRK4Coeffs.beta;
    gamm = p.ETDRK4Coeffs.gamm;
    zeta = p.ETDRK4Coeffs.zeta;
    expMuDt = p.ETDRK4Coeffs.expMuDt;

    % Sub-steps. Some corrections may be required to update
    % time between intermediate steps.

    % 1.
    RHS0 = rightHandSide(p, sol);

    % 2. an's (in Kassam+Trefethen, "n" is time-step).
    an  = expMuDt.*sol + zeta.*RHS0;
    RHS1 = rightHandSide(p, an);

    % 3. bn's and NA2
    bn  = expMuDt.*sol + zeta.*RHS1;
    RHS2 = rightHandSide(p, bn);

    % 4. cn's and NA3
    cn  = expMuDt.*an + zeta.*(2*RHS2-RHS0);
    RHS3 = rightHandSide(p, cn);

    % Time-step
    sol = expMuDt.*sol +   alph .*   RHS0 ...
                       + 2*beta .* ( RHS1 + RHS2 ) ...
                       +   gamm .*   RHS3;

end
