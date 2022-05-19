function [xx,yy,u] = SpectralCollocation_2D_BVP(alpha,beta,gamma,f,N,a,b,BC_type)
% Function to solve the 2D helmholtz equation for ocean acoustic
% propagation given by the form:
% uxx + uyy + a(x,y)ux + b(x,y)uy + c(x,y)u = f(x,y)
% on a square (a,b) x (a,b)
%
% Inputs: 
%        alpha,beta,gamma = variable coefficient function handles
%        f = forcing term function handle
%        N = number of collocation points along one direction
%        a,b = left and right endpoints of x-dimension (respectively)
%        BC_type = type of Boundary Condition
% Outputs:
%        xx,yy = spectral collocation mesh
%        u = approximation to the solution of the
% Authors: Derrick Choi and Ryan Stewart
%
% Note: Method uses Chebyshev spectral collocation
%
%       cheb.m function is from the text Spectral Methods in Matlab by
%       Trefethen. Method for obtaining interior nodes are adapted from
%       poisson_spec code given by Professor Gillman
%       
%       (modification for boundary conditions done by user in separate
%       subfunction and currently assumes one kind of Boundary condition for all
%       sides of the square)


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

% Determine interior vs exterior points
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
Ext_pt = Ext_pt(order_ccw);

% Contributions of variable coefficients
axy = diag(alpha(xx,yy));
bxy = diag(beta(xx,yy));
cxy = diag(gamma(xx,yy));

% check_partial_deriv(xx,yy,du_dx,du_dy,Ext_pt,N,du_dxdx,du_dydy);

% Equations for BCs
switch BC_type
    case 'Dirichlet'    
        [F,fext] = getBC_Eqs(xx,yy,du_dx,du_dy,Ext_pt,N,BC_type);
    case 'Robin'
        [F,fext] = getBC_Eqs(xx,yy,du_dx,du_dy,Ext_pt,N,BC_type);
end

% Equations for interior points
A = du_dxdx+du_dydy+axy*du_dx+bxy*du_dy+cxy;

% Right hand side
fint = f(xx(Int_pt),yy(Int_pt));
b = [fext;fint];

% Solve linear system
FA = [F;A(Int_pt,:)];
u = FA\b;

end

function [F,fext] = getBC_Eqs(xx,yy,dudx,dudy,Ext_pt,N,BC_type)
% Inputs:
%        xx,yy = 2D grid points flattened as long 1D vectors
%        dudx, dudy = partial derivative matrices
%        Ext_pt = exterior point indices
%        N = number of points on a side
%        BC_type = type of BC
% Outputs: 
%        F = Coefficient Matrix for boundary node equations
%        fext = right hand side for boundary nodes
%
Npt = length(xx);

% Matrices
switch BC_type
    case 'Dirichlet'
        F = zeros(length(Ext_pt),Npt);
        F(:,Ext_pt) = eye(length(Ext_pt));
    case 'Neumann'
        % This BC does not work right now
        idx = reshape(1:N,N/4,4);
        D_bottom = dudy(Ext_pt(idx(:,1)),:);
        D_right = dudx(Ext_pt(idx(:,2)),:);
        D_top = dudy(Ext_pt(idx(:,3)),:);
        D_left = dudx(Ext_pt(idx(:,4)),:);
        F = [D_bottom;D_right;D_top;D_left];

    case 'Robin'

        F = zeros(length(Ext_pt),Npt);
        F(:,Ext_pt) = eye(length(Ext_pt));
        
        nodeperside = length(Ext_pt)/4;
        
        % indices for bottom, right, top, and left sides of the square
        b = Ext_pt(1:N-1);
        r = Ext_pt(N:N+nodeperside-1);
        t = Ext_pt(N+nodeperside:N+2*nodeperside-1);
        l = Ext_pt(N+2*nodeperside:end);

        % for the specific problem in the paper, these constants are
        % multiplied to the derivative matrices
        D_bottom = -1/(2*pi)*dudy(b,:);
        D_right = -1/pi*dudx(r,:);
        D_top = -1/(2*pi)*dudy(t,:);
        D_left = -1/pi*dudx(l,:);
        
        % Combine corresponding dirichlet and Neumann entries
        F(1:nodeperside,:) = F(1:nodeperside,:) +D_bottom;
        F(nodeperside+1:2*nodeperside,:) = F(nodeperside+1:2*nodeperside,:)+D_right;
        F(2*nodeperside+1:3*nodeperside,:) = F(2*nodeperside+1:3*nodeperside,:)+D_top;
        F(3*nodeperside+1:end,:) = F(3*nodeperside+1:end,:)+D_left;
        
end

% Values of boundary condition (need to change depending on the problem you
% consider)
%     fext = exp(xx(Ext_pt)+yy(Ext_pt)/2);
    fext = zeros(length(Ext_pt),1);
end

function check_partial_deriv(xx,yy,dudx,dudy,ext,N,dudxdx,dudydy)

% function to check partial derivatives along side of a square
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