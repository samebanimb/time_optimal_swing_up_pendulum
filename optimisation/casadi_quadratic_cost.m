addpath('path/to/casadi')
import casadi.*

n = 4; % state dimension
m = 1; % input dimension 
param.m = 0.127;
param.M = 0.94;
param.l = 0.2066;
param.g = 9.0379;
param.Am= 8.28866833710732;
param.Beq= 1.16342493738728;


x_eq = [0; pi; 0; 0]; % eqilibilium point
u_eq = zeros(m,1);

T  = 1.8;
dt = 0.01;
N  = ceil(T/dt);     % prediction horizon

Q = 10^6*diag([1;1;0;0])+diag([0;0;1;1]);
R = 10^-6;        % input cost

% Initial guess for input
u0 = zeros(m*N,1);
x0 = zeros(n*(N+1),1);
x0(1:n) = [0; 0; 0; 0];
for k=1:N
     x0(n*k+1:n*(k+1)) = dynamic(dt,x0(n*(k-1)+1:n*k),u0(k),param);
end

% constraints
xbound = 0.3;
ubound = 5;

lbx0 = x0(1:n);
lbx  = repmat([-xbound; -Inf; -Inf; -Inf],[N,1]);
lbx(end-n+1:end) = x_eq; % end constraint
lbu  = repmat(-ubound,[N,1]);

ubx0 = x0(1:n);
ubx  = repmat([xbound; Inf; Inf; Inf],[N,1]);
ubx(end-n+1:end) = x_eq;
ubu  = repmat(ubound,[N,1]);

lb   = [lbx0;lbx;lbu];
ub   = [ubx0;ubx;ubu];

con_bound=zeros(N*n,1);
con_lb=con_bound;
con_ub=con_bound; % same bounds enforces equality constraints g(x)=0

% optimization setup
y   = MX.sym('y',N*m+(N+1)*n);
obj = costfunction(N,y,x_eq,u_eq,Q,R,n,m,dt);
con = nonlinearconstraints(N,dt,y,n,m,param);
nlp = struct('x',y,'f',obj,'g',con);
solver = nlpsol('solver','ipopt',nlp);
 
% solve OCP
y_init = [x0;u0];
res = solver('x0' , y_init,...  % solution guess
             'lbx', lb,...      % lower bound on y
             'ubx', ub,...      % upper bound on y
             'lbg', con_lb,...  % lower bound on g(y)
             'ubg', con_ub);    % upper bound on g(y)
         
y_opt = full(res.x); 
x_opt = y_opt(1:n*(N+1));
u_opt = y_opt(n*(N+1)+1:end);
xi    = x_opt(1:4:end);
theta   = x_opt(2:4:end);
dxidt = x_opt(3:4:end);
dthetadt= x_opt(4:4:end);
t=linspace(0,T+dt,N+1);

% plot
tuGreen    =[116,181,0]/255;
tuOrange   =[255,153,0]/255;
tuGray     =[86,86,86]/255;
lightGray  =[170,170,170]/255;
lighterGray=[220,220,220]/255;

%% figures
figure
subplot(3,1,1)
box on
plot(0:dt:T,xi(1:end),'color',tuGreen,'LineWidth',0.8)
yline(xbound,'color',tuGray,'LineStyle','--','LineWidth',0.8)
yline(-xbound,'color',tuGray,'LineStyle','--','LineWidth',0.8)
xlabel('$t$','interpreter','Latex','fontsize',15);
ylabel('$\xi$','interpreter','Latex','fontsize',15);
set(gca,'ylim',xbound*[-1.2,1.2],...
    'TickLabelInterpreter','Latex','fontsize',14)
xlim([0,1.9])

subplot(3,1,2)
box on
plot(0:dt:T,theta(1:end),'color',tuOrange,'LineWidth',0.8)
xlabel('$t$','interpreter','Latex','fontsize',15);
ylabel('$\varphi$','interpreter','Latex','fontsize',15);
set(gca,'ylim',1.2*[-pi,pi],'ytick',-pi:pi/2:pi,'yticklabel',...
    {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'},...
    'TickLabelInterpreter','Latex','fontsize',14)
xlim([0,1.9])

subplot(3,1,3)
box on
stairs(0:dt:T-dt,u_opt,'color','b','LineWidth',0.8)
yline(ubound,'color',tuGray,'LineStyle','--','LineWidth',0.8)
yline(-ubound,'color',tuGray,'LineStyle','--','LineWidth',0.8)
xlabel('$t$','interpreter','Latex','fontsize',15);
ylabel('$u$','interpreter','Latex','fontsize',15);
set(gca,'ylim',ubound*[-1.2,1.2],...
    'TickLabelInterpreter','Latex','fontsize',14)
xlim([0,1.9])

set(gcf,'color','w');
set(gcf,'position',[0 0, 800 500])

%% functions
function xdot = pendulum(x,u,param)
    % System dynamics
    m=param.m;
    M=param.M;
    l=param.l;
    g=param.g;
    Am=param.Am;
    Beq=param.Beq;
    
    xdot = [x(3);
            x(4);
            (m*l*sin(x(2))*x(4)^2+m*g*sin(x(2))*cos(x(2))+(-Am*x(3)+Beq*u))/(M + m*sin(x(2))^2);
            -(m*l*sin(x(2))*cos(x(2))*x(4)^2+(m+M)*g*sin(x(2))+cos(x(2))*(-Am*x(3)+Beq*u))/(l*(M + m*sin(x(2))^2))];
end

function x_new=dynamic(dt,x,u,param)
    % runge kutta forth order
    k1=pendulum(x,u,param);
    k2=pendulum(x+dt/2*k1,u,param);
    k3=pendulum(x+dt/2*k2,u,param);
    k4=pendulum(x+dt*k3,u,param);
    x_new=x+dt/6*(k1+2*k2+2*k3+k4);
end

function cost = costfunction(N,y,x_eq,u_eq,Q,R,n,m,dt)
    % Formulate the cost function to be minimized
    cost = 0;
    x=y(1:n*(N+1));
    u=y(n*(N+1)+1:end);
    
    % Build the cost by summing up the stage cost
    for k=1:N
        x_k=x(n*(k-1)+1:n*k);
        u_k=u(m*(k-1)+1:m*k);
        cost = cost + dt*runningcosts(x_k,u_k,x_eq,u_eq,Q,R);
    end    
end

function cost = runningcosts(x,u,x_eq,u_eq,Q,R)
    % running cost   
    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
end

function con = nonlinearconstraints(N,dt,y,n,m,param)
    % Introduce the nonlinear constraints
    x=y(1:n*(N+1));
    u=y(n*(N+1)+1:end);
    con = [];

    % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);
        u_k=u((k-1)*m+1:k*m);

        ceqnew=x_new - dynamic(dt,x_k,u_k,param);
        con = [con; ceqnew];
    end
end
