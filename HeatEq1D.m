function [tout,xout,u] = HeatEq1D(K,t0,N,a,b,NT,IC,BC)
% function to solve 1D=Heat equation on interval [a,b] fro t = 0 to t = T
% Heat Eq: ut = -k(x)u_xx
%
% Inputs:
% 
% Outputs:
%
% 
%


% Spectral collocation stuff
[D,x] = cheb(N-1);

% map to interval of intereest
x = (b-a)/2*x+(b+a)/2;
J = (b-a)/2; % jacobian
D = 1/J*D; %derivative operator on interval defined by (a,b) 

D2 = D*D; % second derivative on (a,b)

A = diag(K(x(2:end-1)))*D2(2:end-1,2:end-1);

eigval = eig(A);
dt = min(-2./eigval);

tout = t0:dt:dt*(NT-1);

% Solution output initialization
u = zeros(N,length(tout));

% Apply initial condition
u(:,1) = [IC(x)];

% Implicit Euler Time stepping
for i = 2:NT
    u(2:end-1,i) = inv(eye(size(A))-dt*A)*u(2:end-1,i-1);
    u(1,i)= BC(1);
    u(end,i) = BC(2);
end

% RK4 time-stepping scheme
% for i = 2:length(tout)
%     K1 = dt*A*u(2:end-1,i-1);
%     K2 = dt*A*(K1/2+u(2:end-1,i-1));
%     K3 = dt*A*(K2/2+u(2:end-1,i-1));
%     K4 = dt*A*(K3+u(2:end-1,i-1));
% 
%     u(2:end-1,i) = u(2:end-1,i-1)+1/6*(K1+2*K2+2*K3+K4);
%     % Append Dirichlet BCs
%     u(1,i) = BC(1);
%     u(end,i) = BC(2);
% end

xout = x;
end

