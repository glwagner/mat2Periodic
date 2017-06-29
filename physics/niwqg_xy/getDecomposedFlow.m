function [psiqh, psiwh] = getDecomposedFlow(p, sol)

qh = sol(:, :, 2);
Ah = sol(:, :, 1);

Ax = ifft2(1i*p.KK.*Ah);
Ay = ifft2(1i*p.LL.*Ah);
A  = ifft2(Ah);

psiwh = - p.kappa^4./(2*p.f*p.kay2)*( ...
    p.KK.*fft2(conj(A).*Ay) - p.LL.*fft2(conj(A).*Ax)) ...
    - p.kappa^4/(4*p.f)*fft2(conj(A).*A); 

psiqh = -qh ./ p.kay2;
