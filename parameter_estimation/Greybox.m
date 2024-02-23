l =0.3365;
M = 0.94;
m = 0.127;
g = 9.8;
Beq = 1.7235;
Am = 7.72356;
nl_model = idnlgrey('pendulum',[2 1 4], [m,M,l,g,Am,Beq], zeros(4,1),0, 'Name', 'nonlinear_pendulum');

z = iddata([xi_measured(2,:)' theta_measured(2,:)' ],input(2,:)',0.01); 

nl_model.Parameters(1).Fixed = True;
nl_model.Parameters(2).Fixed = True;

opt = nlgreyestOptions;
opt.Display = 'on';

est = nlgreyest(z,nl_model,opt);
est.Name = 'estimated';
%%

function [xdot,y] = pendulum(t,x,u,m,M,l,g,Am,Beq,varargin)
    % System dynamics

    y = [x(1);x(2)];
    xdot = zeros(4,1);
    xdot(1) = x(3);
    xdot(2) = x(4);
    xdot(3)= (m*l*sin(x(2))*x(4)^2+m*g*sin(x(2))*cos(x(2))+...
    (Beq*u-Am*x(3))) /(M + m*sin(x(2))^2);
    xdot(4)= -(m*l*sin(x(2))*cos(x(2))*(x(4)^2)+g*sin(x(2))*(m+M)+...
    (Beq*u-Am*x(3))*cos(x(2)))/(l*(M + m*sin(x(2))^2));

end