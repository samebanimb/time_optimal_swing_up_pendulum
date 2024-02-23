addpath('path/to/casadi')
import casadi.*

% state dimension
n = 4; 

% system parameter
m = 0.127;
M = 0.94;
l = 0.2066;
g =9.03788858603972;
Am= 8.28866833710732;
Beq = 1.16342493738728;
    
% physical constraints
xbound = 0.5;
ubound = 5;



N = 100;
opti = casadi.Opti(); % optimization problem
X = opti.variable(n,N+1); % state trajectory
U = opti.variable(1,N); % control trajectory (throttle)
T = opti.variable(); % final time


opti.minimize(T); % minimize time
f = @(x,u) [x(3);
            x(4);
            (m*l*sin(x(2))*x(4)^2+m*g*sin(x(2))*cos(x(2))+(-Am*x(3)+Beq*u))/(M + m*sin(x(2))^2);
        -(m*l*sin(x(2))*cos(x(2))*x(4)^2+(m+M)*g*sin(x(2))+cos(x(2))*(-Am*x(3)+Beq*u))/(l*(M + m*sin(x(2))^2))]; % dx/dt = f(x,u)
dt = T/N;

for k=1:N % loop over control intervals
k1 = f(X(:,k), U(:,k));
k2 = f(X(:,k)+dt/2*k1, U(:,k));
k3 = f(X(:,k)+dt/2*k2, U(:,k));
k4 = f(X(:,k)+dt*k3, U(:,k));
x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

x = X(1,:);
alpha = X(2,:);
x_dot = X(3,:);
alpha_dot = X(4,:);

%% path constraints
opti.subject_to(-5<=U<=5); % control is limited
opti.subject_to(-0.4<=x<=0.4); % path is limited

%% boundary conditions
opti.subject_to(x(1)==0); % start at position 0 
opti.subject_to(alpha(1)==0); % start from angle 0
opti.subject_to(x_dot(1)==0); % start from stand-still
opti.subject_to(alpha_dot(1)==0); % start at angular speed 0
opti.subject_to(-2.5<=alpha_dot(N+1)<=2.5); % end at angular speed <|2.5| 
opti.subject_to(alpha(N+1)==pi);% finish angle at position pi

%% misc. constraints 
opti.subject_to(T>=0); % Time must be positive

%% initial values for solver 
opti.set_initial(alpha,zeros(1,N+1));
opti.set_initial(alpha_dot,zeros(1,N+1));
opti.set_initial(x,zeros(1,N+1));
opti.set_initial(x_dot, zeros(1,N+1));
opti.set_initial(T,4);
opti.set_initial(U, zeros(1,N));

%% solve NLP 
opti.solver('ipopt'); % set numerical backend
sol = opti.solve(); % actual solve
%%
phi = sol.value(alpha);
dphidt = sol.value(alpha_dot);
xi = sol.value(x);
dxidt = sol.value(x_dot);
T = sol.value(T);
u_opt = sol.value(U);
dt = T/N;
U_opt_3= [(0:dt:T-dt)',u_opt'];
%%
% plot
tuGreen    =[116,181,0]/255;
tuOrange   =[255,153,0]/255;
tuGray     =[86,86,86]/255;
lightGray  =[170,170,170]/255;
lighterGray=[220,220,220]/255;
%%
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
plot(0:dt:T,phi(1:end),'color',tuOrange,'LineWidth',0.8)
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
