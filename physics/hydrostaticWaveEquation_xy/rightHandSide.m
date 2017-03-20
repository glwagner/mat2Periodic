function  RHS = rightHandSide(p, sol)

% -------------------------------------------------------------------------------- 
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

% Get the streamfunction psih from the vorticity qh and set domain mode to zero.
psih = -sol./p.kay2;

% We need dx(psi) and dy(psi). Both are real.
psix = real(ifft2(1i*p.KK.*psih));
psiy = real(ifft2(1i*p.LL.*psih));

% Get vorticity q in real space.  Recall that q is real.
q = real(ifft2(sol));

% Assemble RHS.q. Using J(psi,q) = dy(q psi_x) - dx(q psi_y)
RHS = 1i*p.KK.*fft2(psiy.*q) - 1i*p.LL.*fft2(psix.*q);

% Dealias 
RHS = p.filt.*RHS;
