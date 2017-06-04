function [ daN ] = HIVDelay_RHS( t,aN,kernel )
global na nc c rv ru gam dela delc delu del p r N tau1 tau2 s B0 Bt1 Bt1t2 Q H
global a1 p1 a2 p2
global u11 o11 u12 o12 u21 o21 u22 o22

tj = [0:N]*(r/N);

V = aN(1);
A = aN(2);
C = aN(3);
T = aN(4);
X = A + C + T;

if strcmp(kernel, 'dirac')
    aNt1 = Bt1*aN;
    aNt1t2 = Bt1t2*aN;
    At1 = aNt1(2);
    At1t2 = aNt1t2(2);
elseif strcmp(kernel, 'hat')
    aNt1 = hatKernel(-tau1,s).*aN(2:4:end)';
    At1 = sum(2*r/N * aNt1);
    aNt1t2 = hatKernel(-tau1-tau2,s).*aN(2:4:end)';
    At1t2 = sum(2*r/N * aNt1t2);
elseif strcmp(kernel, 'gamma')
    [lB, uB, gKernel1] = gammaKernel(a1,p1);
    aNt1 = gKernel1.*aN(2:4:end)';
    aNt1 = fliplr(aNt1);
    At1 = trapz(tj(lB):(r/N):tj(uB), aNt1(lB:uB));
    
    [lB, uB, gKernel2] = gammaKernel(a2,p2);
    aNt1t2 = gKernel2.*aN(2:4:end)';
    aNt1t2 = fliplr(aNt1t2);
    At1t2 = trapz(tj(lB):(r/N):tj(uB), aNt1t2(lB:uB));
elseif strcmp(kernel, 'bimod')
    [lB, uB, biKernel1] = biKernel(u11,o11,u12,o12);
    aNt1 = biKernel1.*aN(2:4:end)';
    aNt1 = fliplr(aNt1);
    At1 = trapz(tj(lB):(r/N):tj(uB), aNt1(lB:uB));
    
    [lB, uB, biKernel2] = biKernel(u21,o21,u22,o22);
    aNt1t2 = biKernel2.*aN(2:4:end)';
    aNt1t2 = fliplr(aNt1t2);
    At1t2 = trapz(tj(lB):(r/N):tj(uB), aNt1t2(lB:uB));
end

daN = zeros(size(aN));
daN(1) = -c*V + na*At1 + nc*C - p*V*T;
daN(2) = (rv-dela-del*X)*A - gam*At1t2 + p*V*T;
daN(3) = (rv-delc-del*X)*C + gam*At1t2;
daN(4) = (ru-delu-del*X-p*V)*T;
fourMain = daN(1:4);

daN = Q\(B0'*fourMain + H*aN);
end