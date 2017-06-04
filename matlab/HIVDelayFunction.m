function X = HIVDelayFunction(ti,q)
global na nc c rv ru gam dela delc delu del p S s tau1 tau2 n r N tf B0 Bt1 Bt1t2 Q H 
global a1 p1 a2 p2
%----------- Optimal Invitro Model Parameter Values (Table 4) -------------
kernel = 'gamma';
if length(q) == 3
    p     = q(1);
    tau1  = q(2);
    tau2  = q(3);
elseif length(q) == 2 && ~strcmp(kernel, 'gamma')
    p     = q(1);
    tau1  = q(2);
    tau2  = 0;
elseif length(q) == 1 && ~strcmp(kernel, 'gamma')
    p     = q;
    tau1  = 0;
    tau2  = 0;
elseif length(q) == 1 && strcmp(kernel, 'gamma')
    p     = 8.5e-8;
    tau1  = 24;
    tau2  = 3;
    a1    = q;
    p1    = 29;
    a2    = 1;
    p2    = 3;
elseif length(q) == 2 && strcmp(kernel, 'gamma')
    p     = 8.5e-8;
    tau1  = 24;
    tau2  = 3;
    a1    = q(1);
    p1    = q(2);
    a2    = 1;
    p2    = 3;
elseif length(q) == 4
    p     = 8.5e-8;
    tau1  = 24;
    tau2  = 3;
    a1    = q(1);
    p1    = q(2);
    a2    = q(3);
    p2    = q(4);
end

na    = 0.1194;
nc    = 1.6644e-6;
c     = 0.12;
rv    = 0.035;
ru    = 0.035;
gam   = 8.7625e-4;
dela  = 0.0775;
delc  = 0.025;
delu  = 0.016;
del   = 5.4495e-13;
S     = 0;

s = 4.9;
n = 4;
r = 50;
% if tau1+tau2 == 0
%     r = 24;
% else
%     r = tau1+tau2;
% end
N = 32;
tf = 600;

B0 = B(0);
Bt1 = B(-tau1);
Bt1t2 = B(-tau1-tau2);

% ------------------------- Creating Q matrix -----------------------------
Q = zeros(N+1,N+1);
[m,~] = size(Q);
Q(1:(m+1):end)     = 2/3;  % setting diagonal to 2/3
Q(2:(m+1):end)     = 1/6;  % setting the lower diagonal to 1/6
Q((m+1):(m+1):end) = 1/6;  % setting the upper diagonal to 1/6
Q(1,1)             = N/r + 1/3;  % resetting first element
Q(N+1,N+1)         = 1/3;  % resetting last element
Q = (r/N)*Q;  % multiplying by scaler factor
Q = kron(Q,eye(n));  % Kron multiplication

% --------------- Inner Product of B(0) and (f2(t),0) ---------------------
f2 = [0; 0; 0; S];
h = B0'*f2;
F = Q\h;

% ----------------------- Solving for alpha^N -----------------------------
% Defining eta = [V(0), A(0), T(0), C(0)]' and phi which equals eta at t=0
eta = [0; 1.5e5; 0; 1.35e6];
phi = [0, 1.5e5, 0, 1.35e6];

% Solving for h(eta,phi)
hInt = zeros(n*(N+1),1);
hInt = [hInt(1:4)+(r/(2*N)); hInt(5:end-4)+(r/N); hInt(end-3:end)+(r/(2*N))];
phi = kron(ones(1,N+1),phi);
hInt = hInt.*phi';
hN = B0'*eta + hInt;

% Matrix algebra to solve Q*aN = hN for aN
aN0 = Q\hN;  % coordinate vector to the initial condition

% --------------- Solving h(eta,phi) integral by hand ---------------------
H = zeros(N+1,N+1);
[m,~] = size(H);
H(2:(m+1):end)     = 1/2;  % setting the lower diagonal to 1/2
H((m+1):(m+1):end) = -1/2;  % setting the upper diagonal to 1/2
H(1,1)             = 1/2;  % resetting first element
H(N+1,N+1)         = -1/2;  % resetting last element
H = kron(H,eye(n));  % Kron multiplication

% -------------------- Solve with Runge-Kutta -----------------------------
delt = 0.01;
t = (0:delt:tf)';
nt = numel(t);
aN = zeros(n*(N+1),nt);
aN(:,1) = aN0;
for i = 1:nt-1
    k1 = delt*HIVDelay_RHS(t(i),aN(:,i), kernel);
    k2 = delt*HIVDelay_RHS(t(i)+1/2*delt,aN(:,i)+1/2*k1, kernel);
    k3 = delt*HIVDelay_RHS(t(i)+1/2*delt,aN(:,i)+1/2*k2, kernel);
    k4 = delt*HIVDelay_RHS(t(i)+delt,aN(:,i)+k3, kernel);
    aN(:,i+1) = aN(:,i)+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
end

%------------- Get X at experimental observation times -------------------
for i = 1:length(ti)
    timeInHours = ti(i)*24/delt + 1;
    timeInHours = round(timeInHours);
    V(i) = aN(1,timeInHours);
    A(i) = aN(2,timeInHours);
    C(i) = aN(3,timeInHours);
    T(i) = aN(4,timeInHours);
    X(i) = A(i)+C(i)+T(i);
end

y = [0 1.34e6;3 3.71e6;5 7.53e6;7 1.52e7;10 8.04e6;
     12 9.68e5;14 2.46e6;17 3.26e6;21 1.08e7;24 2.50e7];
 
ti = y(:,1);

% Plotting X(t) with different gamma dists and plotting gamma dists
plots = false;
if kernel == 'gamma' & plots == true;
%     subplot(2,2,1)
%     plot(ti,X,'LineWidth',2)
%     set(gca,'fontsize',16)
%     hold on
%     scatter(ti,X)
%     set(gca,'fontsize',16)
%     xlabel('Days')
%     ylabel('A + C + T')
%     title('Totlal number of cells')


    figure
    scatter(ti, y(:,2))
    set(gca,'fontsize',16)
    hold on
    plot(t/24,aN(2,:)+aN(3,:)+aN(4,:),'LineWidth',2)
    set(gca,'fontsize',16)
    xlabel('Days')
    ylabel('A + C + T')
    title('Total Cells Versus Time')
    
    figure
    subplot(2,2,1)
    plot(t/24,aN(1,:),'LineWidth',2)
    set(gca,'fontsize',16)
    xlabel('Days')
    ylabel('Number of cells')
    title('V')
    
%     subplot(2,2,2)
%     plot(ti,A,'LineWidth',2)
%     set(gca,'fontsize',16)
%     xlabel('Days')
%     ylabel('Number of cells')
%     title('Acutely infected cells')
    
    subplot(2,2,2)
    plot(t/24,aN(2,:),'LineWidth',2)
    set(gca,'fontsize',16)
    xlabel('Days')
    ylabel('Number of cells')
    title('Acutely Infected Cells A')
    
%     subplot(2,2,3)
%     plot(ti,C,'LineWidth',2)
%     set(gca,'fontsize',16)
%     xlabel('Days')
%     ylabel('Number of cells')
%     title('Chronically infected cells')
    
    subplot(2,2,3)
    plot(t/24,aN(3,:),'LineWidth',2)
    set(gca,'fontsize',16)
    xlabel('Days')
    ylabel('Number of cells')
    title('Chronically Infected Cells C')
    
%     subplot(2,2,4)
%     plot(ti,T,'LineWidth',2)
%     set(gca,'fontsize',16)
%     xlabel('Days')
%     ylabel('Number of cells')
%     title('Uninfected target cells')
    
    subplot(2,2,4)
    plot(t/24,aN(4,:),'LineWidth',2)
    set(gca,'fontsize',16)
    xlabel('Days')
    ylabel('Number of cells')
    title(' Uninfected Target Cells T')
end