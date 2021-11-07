
% The drift part of the stochastic system

function dx = mydrift(x,u)
global A B Q R n3 n4 n m;
dx = zeros(m+n+n3+n4,1);

% The origianl system
dx(1:n) = A*x(1:n) + B*u;

% The system of exploration noise
dx(n+1:m+n) = -eye(m)*x(n+1:m+n);

% Compute the z_tilde
Psi = kronv([x(1:n);u;1]);

% Extract $\psi_{t_f}$ in equation (14)
dx(n+m+1:n+m+n4) = sm2vec(Psi*Psi');

% Extract $\xi_{t_f}$ in equation (14)
dx(n+m+n4+1:end) = Psi*(x(1:n)'*Q*x(1:n)+u'*R*u);
end