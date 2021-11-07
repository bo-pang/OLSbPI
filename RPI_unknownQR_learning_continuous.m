
% This scripts learns the near-optimal control policy using the collected
% data

clear;
clc;
global A B C Q R D F n2 n3 n4 n m K p;

 % Read the collected data
load(['data_MIMO.mat']);
[n, m] = size(B);
q1 = length(D);
q2 = length(F);
[~,p] = size(C);
n2 = m+n+1;      
n3 = n2*(n2+1)/2;    
n4 = n3*(n3+1)/2;     

Kest = {K};                % Store the obtained control gain in each iteration
Gest = {};                 % Store the estimated G in each iteration
imax = 20;                 % the total number of learning iterations
i = 1;

% Check stability of the initial control policy
Mat_stab = stochastic_sys_mat(K);
if ~all(eig(Mat_stab)<0)
    disp("Instability");
end

% Find optimal solution using the model-based policy iteration
Ktrue = {K};
Ptrue = {};
for i=1:imax

    % Policy evaluation 
    Ptrue{end+1} = PI(Ktrue{i});

    % Policy improvement
    Ktrue{end+1} = (F{1}'*Ptrue{end}*F{1} + R)\B'*Ptrue{end};
end
Pstar = Ptrue{end};
Kstar = Ktrue{end};

% Learn the near-optimal control policy from data using Algo. 1
inv_Psi_Psi = Psi_Psi\eye(n3);
T_learning = 100;           
for i=1:imax

    % Check the stability of control gain Kest{i}
    Mat_stab = stochastic_sys_mat(Kest{i});
    if  all(all(isnan(Mat_stab))) || ~all(eig(Mat_stab)<0)
        disp("Instability");
        instab = 1;
        break;
    end

    % Solve the ODE on interval [0,T_learning]
    theta0 = sm2vec(zeros(n2));
    [tt,theta] = ode45(@(t,x) Learning_ODE(t,x,Kest{end},inv_Psi_Psi,I_Psi_xx,Psi_r),[0,T_learning],theta0);

    % Extract Ghat, Gest and Kest from theta
    THETA = vec2sm(theta(end,:),n2);
    Ghat = THETA(1:m+n,1:m+n);
    Gest{end+1} = Ghat;
    Kest{end+1} = Ghat(n+1:end,n+1:end)\Ghat(n+1:end,1:n);
end

% Save the learned results into file
save(['results_using_data_MIMO.mat']);

function dx = Learning_ODE(t,x,K,data1,data2,data3)
global n2 m n;
X = vec2sm(x,n2);
dx = data1*(data2*vec(X(1:n,1:n)- X(1:n,n+1:m+n)*K - K'*X(1:n,n+1:m+n)'...
    +K'*X(n+1:m+n,n+1:m+n)*K)+data3);
end
