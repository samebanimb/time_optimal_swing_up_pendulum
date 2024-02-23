T = 3;
   
% initial conditions
ti = [.5;.75;1.25;1.5;1.9];

% state dimnesions
m = length(ti);

dt=0.0;
A=-eye(m)+[zeros(1,m);eye(m-1),zeros(m-1,1)];
b=-dt*ones(m,1);

Aeq = [];
beq = [];

% boxed constraints
lb = zeros(m,1);
ub = 5*ones(m,1);
%% processing
% options for fmincon
tol_opt = 1e-3;
options = optimset('Display','iter','MaxIter',2000,'MaxFunEvals',10000,...
                'TolFun', tol_opt,...
                'Algorithm', 'interior-point',...
                'FinDiffType', 'central');

% Set initial guess
y_init = ti;
objfun = @(y)costfunction(y);
[sol,fval,exitflag,output] = fmincon(objfun,y_init,A,b,Aeq,beq,lb,ub,@(y)nonlcons(y),options);
%% post processing

[X,Tg] = dynamic(sol);
T = sol(end);

%% functions
function xdot = pendulum(t,x,u,m,M,l,g,Am,Beq)
    % System dynamics
    xdot = zeros(4,1);
    xdot(1) = x(3);
    xdot(2) = x(4);
    xdot(3)= (m*l*sin(x(2))*x(4)^2+m*g*sin(x(2))*cos(x(2))+(-Am*x(3)+Beq*u))/(M + m*sin(x(2))^2);
    xdot(4)= -(m*l*sin(x(2))*cos(x(2))*x(4)^2+(m+M)*g*sin(x(2))+cos(x(2))*(-Am*x(3)+Beq*u))/(l*(M + m*sin(x(2))^2));

end

function cost = costfunction(y)
    % Formulate the cost function to be minimized
    [X,~] = dynamic(y);
    cost = y(end);
end

function [c,ceq] = nonlcons(y) 
    [X,Tg] = dynamic(y);
    theta = abs(X(2,end));  
    eps = 0.2618;
    c = [abs(theta-pi)-eps];
    ceq = [];    
end

function [X,Tg] = dynamic(y) 
    %system parameters
    m = 0.127;
    M = 0.94;
    l = 0.2066;
    g =9.03788858603972;
    Am= 8.28866833710732;
    Beq = 1.16342493738728;
    ubound = 5;
    X = zeros(4,1);
    l = length(y);
    Tg = 0;
    for i=1:l
        u = (-1)^i*(ubound);
        [Ti, Xi] = ode45(@(t,x)pendulum (t,x,u,m,M,l,g,Am,Beq), [Tg(end),y(i)], X(:,end));
        Tg = [Tg,Ti(2:end)'];
        X = [X,Xi(2:end,:)'];
    end
end
