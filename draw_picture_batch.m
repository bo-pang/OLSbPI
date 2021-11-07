clc;
clear;
N = 10;
P_error = zeros(N,1);
K_error = zeros(N,1);
J_error = zeros(N,1);
Ptr_error = zeros(N,1);
Ktr_error = zeros(N,1);
Jtr_error = zeros(N,1);
rel_Delta_G = zeros(N,1);
load(['results_using_data_MIMO.mat']);

% % find optimal controller using VI
% PVI = eye(n);
% [tVI,pvi] = ode45(@(t,x) VI_dynamics(t,x),[0,50],vec(PVI));

Jstar = trace(C'*Pstar*C);
for i=1:N
    Ptilde = PI(Kest{i});
    P_error(i) = norm(Ptilde - Pstar,'fro');
    K_error(i) = norm(Kest{i} - Kstar,'fro');
    J_error(i) =  trace(C'*Ptilde*C)-Jstar;
    Ptr_error(i) = norm(Ptrue{i} - Pstar,'fro');
    Ktr_error(i) = norm(Ktrue{i} - Kstar,'fro');
    Jtr_error(i) =  trace(C'*Ptrue{i}*C)-Jstar;
    Gtilde = [Q + A'*Ptilde + Ptilde*A, Ptilde*B;
        B'*Ptilde, R];
    for j=1:q1
        Gtilde(1:n,1:n) = Gtilde(1:n,1:n) + D{j}'*Ptilde*D{j};
    end
    for k=1:q1
        Gtilde(n+1:end,n+1:end) = Gtilde(n+1:end,n+1:end) + F{k}'*Ptilde*F{k};
    end
    rel_Delta_G(i) = norm(Gtilde - Gest{i},'fro')/norm(Gtilde);
    if Jtr_error(i) < 0
        Jtr_error(i) = -Jtr_error(i);
    end
end

subplot(2,2,1);
plot(1:N,K_error,'-bs');hold on;
plot(1:N,Ktr_error,'-rd');
xlabel('Iteration Index');
legend({'$\Vert \hat{K}_i - K^*\Vert_F$','$\Vert K_i - K^*\Vert_F$'},'Interpreter','latex','FontSize',10);
title('(a)');
xlim([1,10]);
xticks([1,3,5,7,10]);

subplot(2,2,2);
plot(1:N,P_error,'-bs');hold on;
plot(1:N,Ptr_error,'-rd');
xlabel('Iteration Index');
legend({'$\Vert \tilde{P}_i - P^*\Vert_F$','$\Vert P_i - P^*\Vert_F$'},'Interpreter','latex','FontSize',10);
title('(b)');
xlim([1,10]);
ylim([0,150]);
xticks([1,3,5,7,10]);

subplot(2,2,3);
plot(1:N,J_error,'-bs');hold on;
plot(1:N,Jtr_error,'-rd');
xlabel('Iteration Index');
legend({'$\Vert \tilde{J}_i - J^*\Vert_F$','$\Vert J_i - J^*\Vert_F$'},'Interpreter','latex','FontSize',10);
title('(c)');
xlim([1,10]);
ylim([0,1.5]);
xticks([1,3,5,7,10]);

subplot(2,2,4);
plot(1:N,rel_Delta_G,'-bs');hold on;
xlabel('Iteration Index');
legend({'$\Vert \Delta{G}_i\Vert_F/\Vert\tilde{G}_i\Vert_F$'},'Interpreter','latex','FontSize',10);
title('(d)');
xlim([1,10]);
ylim([80,110]);
xticks([1,3,5,7,10]);

% function dx = VI_dynamics(t,p)
% global A B D F Q R n q1 q2;
% P = reshape(p,[n,n]);
% DX = Q + A'*P + P*A;
% for i = 1:q1
%     DX = DX + D{i}'*P*D{i};
% end
% Rscr = R;
% for i = 1:q2
%     Rscr = Rscr + F{i}'*P*F{i};
% end
% DX = DX - P*B*(Rscr\B'*P);
% dx = vec(DX);
% end