
% Given a control gian, this function implements the function 
% $\mathcal{A}(\cdot)$ in the paper

function X = stochastic_sys_mat(K)
global A B D F n q1 q2;
X = kron((A-B*K),eye(n)) + kron(eye(n),(A-B*K));
for i=1:q1
    X = X + kron(D{i},D{i});
end
for i=1:q2
    X = X + kron(F{i}*K,F{i}*K);
end
end