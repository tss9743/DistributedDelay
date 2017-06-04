function dsN = Sensitivity_Delay_p1_RHS(t,sN,aN)
global na nc c rv ru gam dela delc delu del p r N Q B0 H
global a1 p1 a2 p2

tj = [0:N]*(r/N);

V = aN(1);
A = aN(2);
C = aN(3);
T = aN(4);

% delay integral
[lB, uB, gKernel1] = gammaKernel(a1,p1);
aNt1 = gKernel1.*aN(2:4:end)';
aNt1 = fliplr(aNt1);
At1 = trapz(tj(lB):(r/N):tj(uB), aNt1(lB:uB));

[lB, uB, gKernel2] = gammaKernel(a2,p2);
aNt1t2 = gKernel2.*aN(2:4:end)';
aNt1t2 = fliplr(aNt1t2);
At1t2 = trapz(tj(lB):(r/N):tj(uB), aNt1t2(lB:uB));

% derivative of gamma wrt alpha
dgKernel = (p1.*a1.^(p1-1).*tj.^(p1-1).*exp(-a1.*tj))./gamma(p1) + ...
           (a1.^(p1).*tj.^(p1-1).*-tj.*exp(-a1.*tj))./gamma(p1);
daNt = dgKernel.*aN(2:4:end)';
daNt = fliplr(daNt);
dAt = trapz(tj(1):(r/N):tj(end), daNt);

% system to solve
dsN = zeros(size(sN));
dsN(1) = (-c - p*T)*sN(1) + nc*sN(3) - p*V*sN(4) + na*At1 + na*dAt;
dsN(2) = p*T*sN(1) + (rv - dela - 2*del*A - del*C - del*T)*sN(2) - del*A*sN(3) + ...
         (-del*A + p*V)*sN(4) - gam*At1t2;
dsN(3) = -del*C*sN(2) + (rv - delc - del*A - 2*del*C - del*T)*sN(3) - del*C*sN(4) + gam*At1t2;
dsN(4) = -p*T*sN(1) - del*T*sN(2) - del*T*sN(3) + (ru - delu - del*A - del*C - 2*del*T - p * V)*sN(4);

% fourMain = dsN(1:4);
% 
% dsN = Q\(B0'*fourMain + H*sN);
end