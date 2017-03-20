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
uh = sol(:, :, 1);
vh = sol(:, :, 2);
hh = sol(:, :, 3);
qh = sol(:, :, 4);

% nonlinear terms, two transforms at a time
u = real(ifft2(uh));
v = real(ifft2(vh));
h = real(ifft2(hh));
q = real(ifft2(qh));

% Linear terms
RHS(:, :, 1) =  p.f0*vh - 1i*p.KK.*hh;
RHS(:, :, 2) = -p.f0*uh - 1i*p.LL.*hh;
RHS(:, :, 3) = -1i*p.cn^2*( p.KK.*uh + p.LL.*vh );

% Get the streamfunction psih from the vorticity qh and set domain mode to zero.
psih = -qh ./ p.kay2;
psix = real(ifft2(1i*p.KK.*psih));
psiy = real(ifft2(1i*p.LL.*psih));

% recall q = psixx+psiyy
psixx = -real(ifft2(p.KK.^2.*psih));
psiyy = -real(ifft2(p.LL.^2.*psih));
psixy = -real(ifft2(p.KK.*p.LL.*psih));

% Nonlinear parts
% u
RHS(:, :, 1) = RHS(:, :, 1) ...
        + fft2(u.*psixy) + fft2(v.*psiyy) ...
		- 1i*p.LL.*fft2(u.*psix) + 1i*p.KK.*fft2(u.*psiy);

% v
RHS(:, :, 2) = RHS(:, :, 2) ...
        - fft2(u.*psixx) - fft2(v.*psixy) ...
		- 1i*p.LL.*fft2(v.*psix) + 1i*p.KK.*fft2(v.*psiy); 
		
% h	
RHS(:, :, 3) = RHS(:, :, 3) ...
		- 1i*p.LL.*fft2(h.*psix) + 1i*p.KK.*fft2(h.*psiy);
		
% q	
RHS(:, :, 4) = -1i*p.LL.*fft2(q.*psix) + 1i*p.KK.*fft2(q.*psiy);

% Dealias 
RHS = p.filt.*RHS;
