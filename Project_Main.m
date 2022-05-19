%% APPM 4660 Project - Main
% Purpose: Code developed for to solve 2D boundary value problems with
% Chebyshev Spectral Collocation
% Author: Derrick Choi
% Date: 4/5/2022
% Last Modified: 4/26/2022

% Things to do: check derivative operations with different kinds of
% functions that exhibit different convergence rates, do a problem with
% variable coefficients, check conditioning of stretching domains, compare
% convergence rates in space and in time
clc;clear;close all;
%% Check Derivative operators

N = 5:5:25;
a = -1;
b = 1;

error = zeros(length(N),4);
error2 = zeros(length(N),4);
error3 = zeros(length(N),4);

g = @(x,y) exp(x).*sin(y);
gx = @(x,y) exp(x).*sin(y);
gy = @(x,y) exp(x).*cos(y);
gxx = @(x,y) exp(x).*sin(y);
gyy = @(x,y) -exp(x).*sin(y);

h = @(x,y) exp(-1./(sin(x/2+y).^2));
hx = @(x,y) (cos(x/2 + y).*exp(-1./sin(x/2 + y).^2))./sin(x/2 + y).^3;
hxx = @(x,y) (cos(x/2 + y).^2.*exp(-1./sin(x/2 + y).^2))./sin(x/2 + y).^6 - (3*cos(x/2 + y).^2.*exp(-1./sin(x/2 + y).^2))./(2*sin(x/2 + y).^4) - exp(-1./sin(x/2 + y).^2)./(2*sin(x/2 + y).^2);
hy = @(x,y) (2*cos(x/2 + y).*exp(-1./sin(x/2 + y).^2))./sin(x/2 + y).^3;
hyy = @(x,y) (4*cos(x/2 + y).^2.*exp(-1./sin(x/2 + y).^2))./sin(x/2 + y).^6 - (6*cos(x/2 + y).^2.*exp(-1./sin(x/2 + y).^2))./sin(x/2 + y).^4 - (2*exp(-1./sin(x/2 + y).^2))./sin(x/2 + y).^2;

r = @(x,y) sin(x.*y+1);
rx = @(x,y) y.*cos(x.*y+1);
ry = @(x,y) x.*cos(x.*y+1);
rxx = @(x,y) -y.^2.*sin(x.*y+1);
ryy = @(x,y) -x.^2.*sin(x.*y+1);

for i = 1:length(N)
    [D,x] = cheb(N(i)-1);
%     x = x(end:-1:1);
    x = (b-a)/2*x+(b+a)/2;
    J = (b-a)/2; % jacobian
    D = 1/J*D;
    D2 = D*D;
    y = x;

    [xx,yy] = meshgrid(x,y);

    xx = xx(:);
    yy = yy(:);

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
    
%     xx = [xx(Ext_pt);xx(Int_pt)];
%     yy = [yy(Ext_pt);yy(Int_pt)];
    % Partial Derivative operators
    du_dydy = kron(eye(N(i)),D2); 
%     du_dxdx = [du_dxdx(Ext_pt,:);du_dxdx(Int_pt,:)];
    du_dxdx = kron(D2,eye(N(i))); 
%     du_dydy = [du_dydy(Ext_pt,:);du_dydy(Int_pt,:)];
    du_dy =   kron(eye(N(i)),D); 
%     du_dy = [du_dy(Ext_pt,:);du_dy(Int_pt,:)];
    du_dx =   kron(D,eye(N(i))); 
%     du_dx = [du_dx(Ext_pt,:);du_dx(Int_pt,:)];

    test_partialx = du_dx*g(xx,yy);
    test_partialy = du_dy*g(xx,yy);
    test_partialxx = du_dxdx*g(xx,yy);
    test_partialyy = du_dydy*g(xx,yy);

    exact_partialx = gx(xx,yy);
    exact_partialy = gy(xx,yy);
    exact_partialxx = gxx(xx,yy);
    exact_partialyy = gyy(xx,yy);

    test_partialx2 = du_dx*h(xx,yy);
    test_partialy2 = du_dy*h(xx,yy);
    test_partialxx2 = du_dxdx*h(xx,yy);
    test_partialyy2 = du_dydy*h(xx,yy);

    exact_partialx2 = hx(xx,yy);
    exact_partialy2 = hy(xx,yy);
    exact_partialxx2 = hxx(xx,yy);
    exact_partialyy2 = hyy(xx,yy);

    test_partialx3 = du_dx*r(xx,yy);
    test_partialy3 = du_dy*r(xx,yy);
    test_partialxx3 = du_dxdx*r(xx,yy);
    test_partialyy3 = du_dydy*r(xx,yy);

    exact_partialx3 = rx(xx,yy);
    exact_partialy3 = ry(xx,yy);
    exact_partialxx3 = rxx(xx,yy);
    exact_partialyy3 = ryy(xx,yy);

    error(i,:) = [ norm(test_partialx-exact_partialx) norm(test_partialy-exact_partialy)...
        norm(exact_partialxx-test_partialxx) norm(exact_partialyy-test_partialyy)];
    error2(i,:) = [ norm(test_partialx2-exact_partialx2) norm(test_partialy2-exact_partialy2)...
        norm(exact_partialxx2-test_partialxx2) norm(exact_partialyy2-test_partialyy2)];
    error3(i,:) = [ norm(test_partialx3-exact_partialx3) norm(test_partialy3-exact_partialy3)...
        norm(exact_partialxx3-test_partialxx3) norm(exact_partialyy3-test_partialyy3)];
end
figure('Position',[200 100 1000 600])
tiledlayout(2,2)
nexttile
loglog(N.^2,error(:,1),'LineWidth',1.5)
hold on
loglog(N.^2,error2(:,1),'LineWidth',1.5)
loglog(N.^2,error3(:,1),'LineWidth',1.5)
xlabel('Number of Nodes','FontSize',12,'Interpreter','latex')
ylabel('$\|u_{exact}-u_{approx}\|_2$','FontSize',12,'Interpreter','latex')
title('${\partial u/\partial x}$','FontSize',14,'Interpreter','latex')
grid on
set(gca,'TickLabelInterpreter','latex')

nexttile
loglog(N.^2,error(:,2),'LineWidth',1.5)
hold on
loglog(N.^2,error2(:,2),'LineWidth',1.5)
loglog(N.^2,error3(:,2),'LineWidth',1.5)
xlabel('Number of Nodes','FontSize',12,'Interpreter','latex')
ylabel('$\|u_{exact}-u_{approx}\|_2$','FontSize',12,'Interpreter','latex')
title('${\partial u/\partial y}$','FontSize',14,'Interpreter','latex')
grid on
set(gca,'TickLabelInterpreter','latex')

nexttile
loglog(N.^2,error(:,3),'LineWidth',1.5)
hold on
loglog(N.^2,error2(:,3),'LineWidth',1.5)
loglog(N.^2,error3(:,3),'LineWidth',1.5)
xlabel('Number of Nodes','FontSize',12,'Interpreter','latex')
ylabel('$\|u_{exact}-u_{approx}\|_2$','FontSize',12,'Interpreter','latex')
title('${\partial^2 u/\partial x^2}$','FontSize',14,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
grid on

nexttile
loglog(N.^2,error(:,4),'LineWidth',1.5)
hold on
loglog(N.^2,error2(:,4),'LineWidth',1.5)
loglog(N.^2,error3(:,4),'LineWidth',1.5)
xlabel('Number of Nodes','FontSize',12,'Interpreter','latex')
ylabel('$\|u_{exact}-u_{approx}\|_2$','FontSize',12,'Interpreter','latex')
title('${\partial^2 u/\partial y^2}$','FontSize',14,'Interpreter','latex')
grid on
set(gca,'TickLabelInterpreter','latex')
sgtitle('\textbf{Convergence Rates for Different Kinds of Functions}','Interpreter','latex','FontSize',16)
lgd = legend('$e^x\sin(y)$','$e^{-1/(\sin(x/2+y)^2)}$','$\sin(xy+1)$','Interpreter','latex','FontSize',12);
lgd.Layout.Tile = 'east';
%% Verify with poisson_spec (it works)
u_exact = @(x,y) exp(x+y/2);
alpha = @(x,y) 0.*x;
beta = @(x,y) 0.*x;
gamma = @(x,y) 0.*x;

a = 0;
b = 1;

f = @(x,y) 1.25*exp(x+y/2);
[xx,yy,u] = SpectralCollocation_2D_BVP(alpha,beta,gamma,f,7,a,b,'Dirichlet');
load('uapptest.mat')

xx = reshape(xx,7,7);
yy = reshape(yy,7,7);
u = reshape(u,7,7);

figure('Name','Test_Poisson')
surfc(xx,yy,u)
xlabel('x','FontSize',12,'Interpreter','latex')
ylabel('y','FontSize',12,'Interpreter','latex')
zlabel('u','FontSize',12,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)
title('\textbf{Solution to $\mathbf{u_{xx}+u_{yy} = 5/4e^{x+y/2}}$}','Interpreter','latex','FontSize',16)
%% Helmholtz Equation for Ocean Acoustics

% Set-up
N = [8 16 22 24 30 32];
a = -1;
b = 1;
alpha = @(x,y) 1*ones(length(x),1);
beta = @(x,y) 1*ones(length(x),1);
gamma = @(x,y) 1*ones(length(x),1);

f = @(x,y) (-5*pi^2 + 1)*sin(pi*x+pi/4).*sin(2*pi*y+pi/4)+...
    pi*cos(pi*x+pi/4).*sin(2*pi*y+pi/4)+2*pi*sin(pi*x+pi/4).*cos(2*pi*y+pi/4);
u_exact = @(x,y) sin(pi*x+pi/4).*sin(2*pi*y+pi/4);

rel_err = zeros(1,length(N));
paper_rel_err = [5.6e-2 2.68e-6 6.24e-11 2.84e-12 8.15e-12 5.51e-12];
for numnodes = 1:length(N)
    % Solve with spectral collocation
    [xgrid,ygrid,u] = SpectralCollocation_2D_BVP(alpha,beta,gamma,f,N(numnodes),a,b,'Robin');
    rel_err(numnodes) = norm(u_exact(xgrid,ygrid)-u)/norm(u_exact(xgrid,ygrid));

end
fprintf(1,'||uex - uapp||/||uex|| = %17.8e\n',norm(u_exact(xgrid,ygrid)-u)/norm(u_exact(xgrid,ygrid)));

% surface plot
xgrid = reshape(xgrid,N(end),N(end));
ygrid = reshape(ygrid,N(end),N(end));
u_grid = reshape(u,N(end),N(end));


% Plot of solution and error
figure('Name','Test_Helmholtz')
surf(xgrid,ygrid,u_grid)
%fsurf(u_exact,[a b a b])

figure('Name','Convergence plot','Position',[400 100 800 600])
semilogy(N.^2,rel_err,'b','LineWidth',1.5);
hold on
semilogy(N.^2,paper_rel_err,'r','LineWidth',1.5)
legend('Relative Error from developed code','Collocation Error from Ma, et al.','Interpreter','latex','FontSize',12)
grid on
xlabel('Number of Nodes','Interpreter','latex','FontSize',12)
ylabel('Relative Error = $\|u_{exact}-u_{approx}\|/\|u_{exact}\|$','FontSize',12,'Interpreter','latex')
title('\textbf{Helmholtz Equation with Robin Boundary Conditions}','Interpreter','latex','FontSize',16)
set(gca,'TickLabelInterpreter','latex','FontSize',14)

%% 1D Heat Equation
close all
x0 = 0;
L = 1;
n = 1:1000;
bk = 320./((2*n-1)*pi);

% set-up
kappa = @(x) 0.0017*ones(length(x),1);
a = 0;
b = 1;
dt = 0.00001;
tend = 10;
IC = @(x) ones(length(x),1);
BC = [0; 0];
N = 10;
NT = 100;

[tout,xout,uout] = HeatEq1D(kappa,0,N,a,b,NT,IC,BC);

%% Animate Figure
test = load('u_test.mat');
f1 = figure('Position',[400 100 800 600]);
%plot(xout,uout(:,1));
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
loops = length(uout(1,:));
M(loops) = struct('cdata',[],'colormap',[]);
f1.Visible = 'off';
mov = zeros(612,775,1,loops,'uint8');
for j = 1:loops
    plot(xout,uout(:,j),'LineWidth',1.5)
    hold on
    plot(test.x,test.u(:,j),'LineWidth',1.5)
    legend('Spectral Collocation','FEM','Interpreter','latex')
    xlabel('x-axis','FontSize',14,'Interpreter','latex')
    ylabel('y-axis','FontSize',14,'Interpreter','latex')
    title('\textbf{Solution to 1D Heat Equation}','Interpreter','latex','FontSize',16)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    drawnow
    M(j) = getframe;
    hold off
    if j == 1
        [mov(:,:,1,j), map] = rgb2ind(M(j).cdata, 256, 'nodither');
    else
        mov(:,:,1,j) = rgb2ind(M(j).cdata, map, 'nodither');
    end
end
f1.Visible = 'on';
movie(M);
% Create gif
imwrite(mov,map,'HeatEQ1D.gif','DelayTime',0,'LoopCount',Inf)
% approx exact
% u_exact = @(x,t) 100-sum(bk.*sin((2*n-1)*pi.*x)/L.*exp(-((2*n-1)*pi/L).^2*0.0017.*t));
% u_ap = zeros(length(xout),length(xout));
% Evaluate Exact
% for j = 1:length(xout)
%     for i = 1:length(tout)
%         u_ap(j,i) = u_exact(xout(j),tout(i));
%     end
% end

%% 2D-Heat Equation
K = @(x,y) ones(length(x),1);
% f = @(x,y,t) 0.*x; % test problem with no forcing
f = @(x,y,t) exp(-t.^2).*sin(pi.*x.*y).*(-2*t+pi^2*(x.^2+y.^2));
a = -1;
b = 1;
N = 5:5:30;

dt = min(1./(0.048*2*N.^4));
t0 = 0;
tend = 2;
% Time Error
t_err = cell(2,length(dt));
con = zeros(length(dt),1);

for m = 1:length(N)
    fprintf('Solving with Spectral Collocation: dt = %f\n',dt)
    [tout,xout,yout,uout,con(m)] = HeatEq2D(K,f,dt,t0,tend,N(m),a,b,'Dirichlet');

    [xgrid,ygrid] = meshgrid(xout,yout);
    % Compare Error
    % u_exact = @(x,y,t) exp(-2*t).*sin(x).*sin(y); % exact solution for test problem
    u_exact = @(x,y,t) exp(-t.^2).*sin(pi*x.*y);
    fprintf('Computing Exact Solution!\n')
    exact = zeros(N(m),N(m),length(tout));
    for k = 1:length(tout)
        for i = 1:length(xout)
            for j = 1:length(yout)
                exact(i,j,k) = u_exact(xout(i),yout(j),tout(k));
            end
        end
    end
    fprintf('Computing Error in Time!\n')
    for id_t = 1:length(tout)
        approx = reshape(uout(:,id_t),N(m),N(m));
        err = reshape(exact(:,:,id_t)-approx,N(m)^2,1);
        t_err{1,m}(id_t,1) = norm(err,2);
    end
    t_err{2,m} = tout;
end

figure('Position',[400 100 1000 600])
tiledlayout(1,2)
nexttile

surf(xgrid,ygrid,exact(:,:,end))
xlabel('x-axis','Interpreter','latex','FontSize',12)
ylabel('y-axis','Interpreter','latex','FontSize',12)
zlabel('u','FontSize',12,'Interpreter','latex','Rotation',0);
set(gca,'TickLabelInterpreter','latex','FontSize',12)
title('\textbf{Exact Solution}','Interpreter','latex','FontSize',14)

nexttile
surf(xgrid,ygrid,reshape(uout(:,end),N(m),N(m)))
xlabel('x-axis','Interpreter','latex','FontSize',12)
ylabel('y-axis','Interpreter','latex','FontSize',12)
zlabel('u','FontSize',12,'Interpreter','latex','Rotation',0);
set(gca,'TickLabelInterpreter','latex','FontSize',12)
title('\textbf{Approximate Solution}','Interpreter','latex','FontSize',14)

sgtitle('\textbf{Comparison of Solutions}','FontSize',16,'Interpreter','latex')
cbar = colorbar;
cbar.Layout.Tile = 'east';
cbar.Label.String = 'u';
cbar.Label.Interpreter = 'latex';
cbar.TickLabelInterpreter = 'latex';
cbar.FontSize = 12;
cbar.Label.Rotation = 0;
% black = [0 0 0];
% gold = [207 184 124]/255;
% colors_p = [linspace(black(1),gold(1),1000)', linspace(black(2),gold(2),1000)', linspace(black(3),gold(3),1000)'];
% colormap(colors_p)
%{
%%
figure('Name','HeatEq2D_TimeError','Position',[400 100 800 600])
for m = 1:size(t_err,2)
semilogy(t_err{2,m},t_err{1,m},'LineWidth',1.5)
hold on
end
xlabel('Time, t','Interpreter','latex','FontSize',14)
ylabel('$\|u_{exact}(x,y)-u_{approx}(x,y)\|_2$','FontSize',14,'Interpreter','latex')
title('$\mathbf{L_2}$\textbf{-error at each time step}, $\Delta \mathbf{t = \frac{10}{N^{4}}}$','Interpreter','latex','FontSize',16)
grid on
% legend('$N = 25$','$N = 100$','$N = 225$','$N = 400$','$N = 625$','$N = 900$','Interpreter','latex','FontSize',12)
%legend('$\Delta t = 1.9753 \cdot 10^{-5}$','$\Delta t = 9.8765 \cdot 10^{-6}$','$\Delta t = 4.9383 \cdot 10^{-6}$'...
%   ,'Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
%%
%}
%%
f2 = figure('Position',[400 100 800 600]);
axis tight manual
ax = gca;
% ax.View = [-12.2,51.22];
% ax.CameraPosition = [-2.311012994955346,-10.599436255823809,14.842021805308326];
% ax.CameraViewAngle = [-10.2541];
zlim([min(uout,[],'all')-0.1 max(uout,[],'all')+0.1])
grid on
ax.NextPlot = 'replaceChildren';
loops = length(uout(1,:));
M(loops) = struct('cdata',[],'colormap',[]);
f2.Visible = 'off';

it = 1:100:loops;
v = VideoWriter('HeatEq2DAnim.mp4');
open(v)
for j = 1:length(it)
    approx_u = reshape(uout(:,it(j)),N,N);
    surf(xgrid,ygrid,approx_u,'EdgeColor','flat')
    xlabel('x-axis','FontSize',14,'Interpreter','latex')
    ylabel('y-axis','FontSize',14,'Interpreter','latex')
    zlabel('u(x,y)','FontSize',14,'Interpreter','latex')
    title('\textbf{Solution to 2D Heat Equation with Forcing}','Interpreter','latex','FontSize',16)
    ax.TickLabelInterpreter = 'latex';
    xlim([a b])
    ylim([a b])
    grid on
    colorbar
    drawnow
    M(j) = getframe;
    im{j} = frame2im(M(j));
    writeVideo(v,M(j))
end
close(v)

for idx = 1:length(it)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,'HeatEQ2D.gif','gif','LoopCount',Inf,'DelayTime',0)
    else
        imwrite(A,map,'HeatEQ2D.gif','gif','WriteMode','append','DelayTime',0)
    end
end