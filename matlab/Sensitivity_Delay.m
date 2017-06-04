% ------------ Get the solutions to the delay system ----------------------
run('HIVDelay.m')

% -------------------- Solve with Runge-Kutta -----------------------------
delt = 0.01;
t = (0:delt:600)';
nt = numel(t);
sN = zeros(4,nt);
for i = 1:nt-1
    k1 = delt*Sensitivity_Delay_p1_RHS(t(i),sN(:,i), aN(:,i));
    k2 = delt*Sensitivity_Delay_p1_RHS(t(i)+1/2*delt,sN(:,i)+1/2*k1, aN(:,i));
    k3 = delt*Sensitivity_Delay_p1_RHS(t(i)+1/2*delt,sN(:,i)+1/2*k2, aN(:,i));
    k4 = delt*Sensitivity_Delay_p1_RHS(t(i)+delt,sN(:,i)+k3, aN(:,i));
    sN(:,i+1) = sN(:,i)+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
end