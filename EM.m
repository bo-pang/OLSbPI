

% Euler-Maruyama method for the stochastic integration
% from ta to tb with time-step dt, initial state x0, dynamics f and g


function [xaugend,u_end,I_Psi_xx,Psi_Psi,Psi_r] = EM(f,g,x0,dt,ta,tb)
global m n n3 n4 K p q1 q2 Noise_Mag;
N = floor((tb-ta)/dt);
x = x0;
u = [];

% Calculate the stochastic integral
for i=1:N
    rv = normrnd(0,1,[p+m+q1+q2,1]);
    u = [u, -K*x(1:n,end) + Noise_Mag*x(n+1:m+n,end)];
    x_next = x(:,end) + f(x(:,end),u(:,end))*dt + ...
        g(x(:,end),u(:,end))*rv*dt^0.5;
    x = [x x_next];
end
result = x;
[~,tN] = size(result);

% Calculate $\zeta_{t_f}$ using the definition of Riemann integral
I_Psi_xx = zeros([n3,n^2]);
for i=1:tN-1
    xx = kron(result(1:n,i),result(1:n,i));
    xx_next = kron(result(1:n,i+1),result(1:n,i+1));
    Psi = kronv([result(1:n,i);u(:,i);1]);
    delta = Psi*(xx_next-xx)';
    I_Psi_xx = I_Psi_xx + delta;
end

Psi_Psi = vec2sm(result(n+m+1:n+m+n4,end),n3);
Psi_r = result(n+m+n4+1:end,end);
xaugend = result(1:m+n,end);
u_end = u(:,end);

end