function [u, v] = getVelocities(p, Ah)

    uh = 1/(p.alpha*p.f0)*( p.sigma*p.KK.*Ah + 1i*p.f0*p.LL.*Ah );
    vh = 1/(p.alpha*p.f0)*( p.sigma*p.LL.*Ah - 1i*p.f0*p.KK.*Ah );
                
    u = ifft2(uh);
    u = real(u + conj(u));

    v = ifft2(vh);
    v = real(v + conj(v));

end
