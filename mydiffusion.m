
% The drift part of the stochastic system

function DD = mydiffusion(x,u)
global C D F n3 n4 m p n q1 q2;

% Compute the gain matrix of the multiplicative noise
dependent = [];
for i=1:q1
    dependent = [dependent, D{i}*x(1:n)];
end
for i=1:q2
    dependent = [dependent, F{i}*u];
end

% DD(1:n,:): the gain matrix of the additive and multiplicative noise
% DD(n+1:n+m): the gain matrix of the exploration noise
% DD(n+m+1:end): the gain matrix of the integral in equation (14)
DD = [C zeros(n,m) dependent;
    zeros(m,p) eye(m) zeros(m,q1+q2); 
    zeros(n3+n4,p+m+q1+q2)];
end