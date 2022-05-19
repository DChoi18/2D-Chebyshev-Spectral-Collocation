function [tout,xout,yout,u,con] = HeatEq2D(K,f,dt,t0,tend,N,a,b,BC_type)
% function to solve 2D-Heat equation on square [a,b] x [a,b] from t = 0 to t ~ T
% Heat Eq: ut - k(x)(u_xx+uyy) = f(x,y)
%
%
% Inputs:
%        K = material thermal diffusivity
%        f = forcing function
%        t0 = start time
%        tend = stop time
%        dt = time step size
%        N = Number of nodes for spatial discretization in one dimension
%        a = "left" end point of side of square
%        b = "left" end point of side of square
%        BC_type = type of Boundary condition on the boundary of domain
% Outputs:
%        tout = vector of times solution is approximated at
%        xout,yout = spatial discretization mesh points
%        u = solution at (x,y) at time = t
%        con = condition number of the matrix to be solved
%
% Note: initial conditions and BCs applied in a subroutine
%
%       cheb.m function is from the text Spectral Methods in Matlab by
%       Trefethen. Method for obtaining interior nodes are adapted from
%       poisson_spec code given 

% Get chebyshev differentiation matrix for interior points
[D,x] = cheb(N-1);
x = x(end:-1:1);
% map to interval of interest
x = (b-a)/2*x+(b+a)/2;
J = (b-a)/2; % jacobian
D = 1/J*D; %derivative operator on interval defined by (a,b) 

% second derivative differentiation matrix
D2 = D*D; % 1D second derivative matrix

% Create tensor product grid
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:);
yy = yy(:);

% Derivative approximations of the left hand side
du_dydy = kron(eye(N),D2);
du_dxdx = kron(D2,eye(N));
du_dy =   -kron(eye(N),D);
du_dx =   -kron(D,eye(N));

% Determine interior vs boundary points
logic_idx = xx==a|xx==b|yy==a|yy==b;
Ext_pt = find(logic_idx);
Int_pt = find(~logic_idx);

% Order exterior points so that it goes from bottom left corner and moves
% counterclockwise
xcenter = (b+a)/2;
ycenter = (b+a)/2;

% reference angle to determine point locations relative to center
theta0 = atan2(yy(1)-ycenter,xx(1)-xcenter);
theta = rem(4*pi+1e-12-theta0+atan2(yy(Ext_pt)-ycenter,xx(Ext_pt)-xcenter),2*pi);
[~,order_ccw] = sort(theta);
Ext_pt = Ext_pt(order_ccw); % gives indexing so that points are ordered from bottom left and moves ccw

% Contributions of variable coefficients
kappa = diag(K(xx,yy));

% Spatial Discretization operator
A = du_dx*kappa*du_dx+du_dy*kappa*du_dy;
Aint = A(Int_pt,:); % only for the interior points

F = zeros(length(Ext_pt),N^2);
F(:,Ext_pt) = eye(length(Ext_pt));

%eigval = eig([F;Aint]);

% Flow: initialize t = 0, advance in time the interior nodes, apply BCs to 
% boundary at each timestep
tout = t0:dt:tend;

u = zeros(length(xx),length(tout));
con = cond([F;Aint]);

u(:,1) = applyICs(xx,yy,Ext_pt,Int_pt,N);

% check_partial_deriv(xx,yy,du_dx,du_dy,Ext_pt,N,du_dxdx,du_dydy);

% Begin time stepping scheme (Euler)
for it = 2:length(tout)
    u(Int_pt,it) = u(Int_pt,it-1) + dt*(Aint*u(:,it-1)+f(xx(Int_pt),yy(Int_pt),tout(it-1)));
    % apply BCs to boundary nodes
    [~,u(Ext_pt,it)] = applyBCs(tout(it),xx,yy,du_dx,du_dy,Ext_pt,N,BC_type);
end

% output grid
xout = x;
yout = y;
end

%% Initial Conditions
function u0 = applyICs(xx,yy,boundary,interior,N)

% function to apply initial conditions 

% IC = @(x,y) sin(y).*sin(x); % initial condition function for test problem
IC = @(x,y) sin(pi*x.*y);
nodeperside = length(boundary)/4;

% identify bottom, right, top, and left indices of square
b = boundary(1:N-1);
r = boundary(N:N+nodeperside-1);
t = boundary(N+nodeperside:N+2*nodeperside-1);
l = boundary(N+2*nodeperside:end);

% apply IC to nodes
u0(1:nodeperside,:) = IC(xx(b),yy(b));
u0(nodeperside+1:2*nodeperside,:) = IC(xx(r),yy(r));
u0(2*nodeperside+1:3*nodeperside,:) = IC(xx(t),yy(t));
u0(3*nodeperside+1:length(boundary),:) = IC(xx(l),yy(l));
u0(length(boundary)+1:length(boundary)+length(interior),:) = IC(xx(interior),yy(interior));

u0 = zeros(length(boundary)+length(interior),1);
u0(boundary) = IC(xx(boundary),yy(boundary));
u0(interior) = IC(xx(interior),yy(interior));
end

%% Boundary Conditions
function [F,fext] = applyBCs(time,xx,yy,dudx,dudy,Ext_pt,N,BC_type)
% Inputs:
%        xx,yy = 2D grid points flattened as long 1D vectors
%        dudx, dudy = partial derivative matrices
%        Ext_pt = exterior point indices
%        N = number of points on a side
%        BC_type = type of BC

Npt = length(xx);

% Matrices
switch BC_type
    case 'Dirichlet'

        F = zeros(length(Ext_pt),Npt);
        F(:,Ext_pt) = eye(length(Ext_pt));
        
        nodeperside = length(Ext_pt)/4;
        
        % indices for bottom, right, top, and left sides of the square
        b = Ext_pt(1:N-1);
        r = Ext_pt(N:N+nodeperside-1);
        t = Ext_pt(N+nodeperside:N+2*nodeperside-1);
        l = Ext_pt(N+2*nodeperside:end);
        
end

% Values of boundary condition (need to change depending on the problem you
% consider)
% fbot = zeros(length(b),1);
% fright = zeros(length(r),1);
% ftop = zeros(length(t),1);
% fleft = zeros(length(l),1);
fbot = -exp(-time^2)*sin(pi*xx(b));
fleft = -exp(-time^2)*sin(pi*yy(l));
ftop = exp(-time^2)*sin(pi*xx(t));
fright = exp(-time^2)*sin(pi*yy(r));

fext = [fbot;fright;ftop;fleft];
end

function check_partial_deriv(xx,yy,dudx,dudy,ext,N,dudxdx,dudydy)

s = ext(1:N);
e = ext(N+1:N+5);
n = ext(N+6:N+11);
w = ext(N+12:end);

u = @(x,y) exp(x).*sin(y);
ux = @(x,y) exp(x).*sin(y);
uy = @(x,y) exp(x).*cos(y);

Dsx = dudx(s,:);Dex = dudx(e,:);
Dnx = dudx(n,:);Dwx = dudx(w,:);

Dsy = dudy(s,:);Dey = dudy(e,:);
Dny = dudy(n,:);Dwy = dudy(w,:);

ff = u(xx,yy);

test{:,1} = Dsx*ff; exact{:,1} = ux(xx(s),yy(s));
test{:,2} = Dwx*ff; exact{:,2} = ux(xx(w),yy(w));
test{:,3} = Dnx*ff; exact{:,3} = ux(xx(n),yy(n));
test{:,4} = Dex*ff; exact{:,4} = ux(xx(e),yy(e));
test{:,5} = Dsy*ff; exact{:,5} = uy(xx(s),yy(s));
test{:,6} = Dwy*ff; exact{:,6} = uy(xx(w),yy(w));
test{:,7} = Dny*ff; exact{:,7} = uy(xx(n),yy(n));
test{:,8} = Dey*ff; exact{:,8} = uy(xx(e),yy(e));
 
for i = 1:length(test)
    err(i) = norm(exact{:,i}-test{:,i});
end

figure
plot(err)
% figure
% plot(ux_ex)
% hold on
% plot(test)
% err = norm(ux_ex-test);
disp(err)
end
