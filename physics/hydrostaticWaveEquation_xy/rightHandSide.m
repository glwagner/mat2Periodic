function  RHS = rightHandSide(p, sol)

% Inputs:
%	p		: parameters
%	sol 	: solution
%
% Output:
%   RHS      : Nonlinear RHS of the PDE being solved.
% 
% For calculating Jacobian terms we use the identities
%
%					 J(a,b)	= a_x b_y - a_y b_x ,
%							= dx(a b_y) - dy(a b_x) , 
%							= dy(b a_x) - dx(b a_y) .
%
% The transform of the Jacobian can therefore be written
%
%					{J(a,b)} = ik{a b_y} - il{a b_x} , 
%							 = il{b a_x} - ik{b a_y} ,
%
% where {a b} denotes the transform of a*b.
%

% Key to the solution
Ah = sol(:, :, 1);
qh = sol(:, :, 2);

% Get the streamfunction psih from the vorticity qh and set domain mode to zero.
psih = -qh ./ p.kay2;

U = -real(ifft2(1i*p.LL.*psih));
V =  real(ifft2(1i*p.KK.*psih));
q = real(ifft2(qh));

Ax = ifft2(1i*p.KK.*Ah);
Ay = ifft2(1i*p.LL.*Ah);
EA = -p.alpha/2 * ifft2(( p.KK.^2+p.LL.^2 + p.kappa^2*(4+3*p.alpha) ).*Ah );
        
Axx = -ifft2(p.KK.^2.*Ah);
Ayy = -ifft2(p.LL.^2.*Ah);
Axy = -ifft2(p.KK.*p.LL.*Ah);

% Construct RHS.
RHS = zeros(p.ny, p.nx, p.nVars);

% For A in parts:
% 1. Advection
RHS(:, :, 1) = RHS(:, :, 1) ...
    -p.invE.*( 1i*p.KK.*fft2(U.*EA) + 1i*p.LL.*fft2(V.*EA) );

% 2. Refraction
RHS(:, :, 1) = RHS(:, :, 1) ...
    -p.invE/p.f0.*(  1i*p.KK.*fft2( (1i*p.sigma*Ax - p.f0*Ay).*q ) ...
                   + 1i*p.LL.*fft2( (1i*p.sigma*Ay + p.f0*Ax).*q ) );

% 3. Middling Jacobian terms
RHS(:, :, 1) = RHS(:, :, 1) + 2i*p.invE*p.sigma/p.f0^2 .* ( ... 
    1i*p.KK.*fft2(  V.*( 1i*p.sigma*Axy - p.f0*Ayy ) ...
                  - U.*( 1i*p.sigma*Ayy + p.f0*Axy ) ) ...
  - 1i*p.LL.*fft2(  V.*( 1i*p.sigma*Axx - p.f0*Axy ) ...
                  - U.*( 1i*p.sigma*Axy + p.f0*Axx ) ) ); 

% q	
RHS(:, :, 2) = -1i*p.KK.*fft2(U.*q) - 1i*p.LL.*fft2(V.*q);

% Dealias 
RHS = p.filt.*RHS;
