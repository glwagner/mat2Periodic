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

Ax = ifft2(1i*p.KK.*Ah);
Ay = ifft2(1i*p.LL.*Ah);
A  = ifft2(Ah);

% Mean and wave-induced mean flow
psih = -qh./p.kay2 ...
    - p.kappa^4./(2*p.f*p.kay2)*( ...
    p.KK.*fft2(conj(A).*Ay) - p.LL.*fft2(conj(A).*Ax)) ...
    - p.kappa^4/(4*p.f)*fft2(conj(A).*A); 

psih(1, 1) = 0;

U = -real(ifft2(1i*p.LL.*psih));
V =  real(ifft2(1i*p.KK.*psih));
Z = -real(ifft2((p.KK.^2+p.LL.^2).*psih));
q = real(ifft2(qh));
        
% Construct RHS.
RHS = zeros(p.ny, p.nx, p.nVars);

% For A:
RHS(:, :, 1) = -p.invL.*(1i*p.KK.*fft2(U.*A) + 1i*p.LL.*fft2(V.*A) ...
    + 1i/2*fft2(A.*Z));
    
% q	
RHS(:, :, 2) = -1i*p.KK.*fft2(U.*q) - 1i*p.LL.*fft2(V.*q);

% Dealias 
RHS = p.filt.*RHS;
