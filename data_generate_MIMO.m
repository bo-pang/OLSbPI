
% This script generates and collects the data to be used in the learing
% phase from the stochastic system

clear;
clc;
global A B D F C Q R n2 n3 n4 n m K p q1 q2 c3 Noise_Mag;
Noise_Mag = 1e2;
% Fix the random seed
rng(1);

% The model of the triple inverted pendulum
A = [zeros(3), eye(3);
    12.54 -8.26 -0.39 -0.043 2.75 -0.36;
    -4.38 36.95 -3.00 0.086 -9.57 2.29;
    -6.82 -22.94 11.93 -0.034 6.82 -2.86;];
B = [zeros(3,2);
    -50.0 6.12;
    174.4 -38.93;
    -124.2 48.62;];

% Set the multiplicative noise
D = {};
F = {};
sigma = 0.01;
D{1} = diag([0,0,0,0,0,sigma]);
F{1} = [zeros(3,2);
    sigma 0;
    zeros(2);];

% Get the dimension of different variables
q1 = length(D);
q2 = length(F);
[n, m] = size(B);

% Set the additive noise
c3 = 0.1;
C = c3*eye(n);
[~,p] = size(C);
Q = eye(n);
R = eye(m);

% Set the initial controller, which is assumed known in the paper but is
% not necessarily optimal.
K = place(A,B,[-6:1:-1]);

n2 = m+n+1;            % Dimension of z in equation (14)
n3 = n2*(n2+1)/2;      % Dimension of $\xi_{t_f}$ in equation (14)
n4 = n3*(n3+1)/2;      % Dimension of $\psi_{t_f}$ in equation (14)

kmax = 1e4;        % Numer of simulation iterations in data collection
k = 1;             
t = 10;        
dt = 0.0001;       % The time-step in the Eular-Maruyama method 
tstep = 0.05;      
% store the state and input trajectory
x_his = [];        
u_his = [];

% Collect data
while k<kmax
    k/kmax         % See the progress of data collection
    if k==1
        xaug0 = [0.1,0.2,0.3,0.4,0.5,0.6,0,0]';
        y0 = zeros(n3+n4,1);
        % Initial state of the integration
        % X0[1:n]: the state of the original system equation (1)
        % X0[n+1:n+m]: the state of the exploration noise equation (10)
        % X0[n+m+1:n+m+n4]: the vectorization of $\psi_{t_f}$ in equation (14)
        % X0[n+m+n4+1:end]: $\xi_{t_f}$ in equation (14)
        X0 = [xaug0;y0];
        % Implement the stochastic integration
        % I_Psi_xx: $\zeta_{t_f}$ in equation (14)
        % Psi_Psi: $\psi_{t_f}$ in equation (14)
        % Psi_r: $\xi_{t_f}$ in equation (14)
        [xaug, u_end, I_Psi_xx, Psi_Psi, Psi_r] = EM(@mydrift,@mydiffusion,X0,dt,0,t);
    else
        t = t + tstep;
        [xaug, u_end, d_I_Psi_xx, d_Psi_Psi, d_Psi_r] = EM(@mydrift,@mydiffusion,X,dt,t-tstep,t);
        I_Psi_xx = I_Psi_xx + d_I_Psi_xx;
        Psi_Psi = Psi_Psi + d_Psi_Psi;
        Psi_r = Psi_r + d_Psi_r;
    end

    X = [xaug;y0];
    k = k + 1;
    x_his = [x_his,xaug(1:n)];
    u_his = [u_his,u_end];


%     % Check the rank condition empirically
%     rk = rank(Psi_Psi/t);
%     if rk ~= n3
%         disp('Rank defficient');
%         break;
%     end
end
% Save the generated data for the usage in the learning phase
save('data_MIMO.mat');
