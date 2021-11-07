clc;
clear;
load("results_using_data_MIMO.mat");

dt = 0.001;
N = 1e4;
x0 = randn(n,1)*10;
dx = zeros(n,1);
rng(1);

% check stability
Mat_stab = stochastic_sys_mat(Kest{1});
if ~all(eig(Mat_stab)<0)
    disp("Instability");
end
x1 = x0;
rv1 = randn(q1,N);
rv2 = randn(q2,N);
rv3 = randn(n,N);
for i=1:N-1
%     i/N
    u = -Kest{1}*x1(:,end);
    dx = (A*x1(1:n,end)+B*u)*dt + C*rv3(:,i)*sqrt(dt)...
        + D{1}*x1(1:n,end)*rv1(:,i)*sqrt(dt) + F{1}*u*rv2(:,i)*sqrt(dt);
    x1 = [x1, x1(:,end)+dx];
end


% check stability
Mat_stab = stochastic_sys_mat(Kest{10});
if ~all(eig(Mat_stab)<0)
    disp("Instability");
end
x2 = x0;
rv1_n = randn(q1,N);
rv2_n = randn(q2,N);
rv3_n = randn(n,N);
for i=1:N-1
%     i/N
    u = -Kest{10}*x2(:,end);
    dx = (A*x2(1:n,end)+B*u)*dt + C*rv3_n(:,i)*sqrt(dt)...
        + D{1}*x2(1:n,end)*rv1_n(:,i)*sqrt(dt) + F{1}*u*rv2_n(:,i)*sqrt(dt);
    x2 = [x2, x2(:,end)+dx];
end

% Plot the trajectory in data collection phase
[~,simN] = size(x_his);

figure(1);
subplot(3,2,1);
plot(tstep*[0:simN-1],x_his(1,:));
xlim([0,500]);
title('$x_1$','Interpreter','latex');
xlabel('Time(s)');
subplot(3,2,2);
plot(tstep*[0:simN-1],x_his(2,:));
xlim([0,500]);
title('$x_2$','Interpreter','latex');
xlabel('Time(s)');
subplot(3,2,3);
plot(tstep*[0:simN-1],x_his(3,:));
xlim([0,500]);
title('$x_3$','Interpreter','latex');
xlabel('Time(s)');
subplot(3,2,4);
plot(tstep*[0:simN-1],x_his(4,:));
xlim([0,500]);
title('$x_4$','Interpreter','latex');
xlabel('Time(s)');
subplot(3,2,5);
plot(tstep*[0:simN-1],x_his(5,:));
xlim([0,500]);
title('$x_5$','Interpreter','latex');
xlabel('Time(s)');
subplot(3,2,6);
plot(tstep*[0:simN-1],x_his(6,:));
xlim([0,500]);
title('$x_6$','Interpreter','latex');
xlabel('Time(s)');

% Plot the trajectory after learning
figure(2)

subplot(3,2,1);
plot(dt*[0:N-1],x1(1,:),'b--');hold on;
plot(dt*[0:N-1],x2(1,:),'r-');
title('$x_1$','Interpreter','latex');
xlabel('Time(s)');

subplot(3,2,2);
plot(dt*[0:N-1],x1(2,:),'b--');hold on;
plot(dt*[0:N-1],x2(2,:),'r-');
title('$x_2$','Interpreter','latex');
xlabel('Time(s)');

subplot(3,2,3);
plot(dt*[0:N-1],x1(3,:),'b--');hold on;
plot(dt*[0:N-1],x2(3,:),'r-');
title('$x_3$','Interpreter','latex');
xlabel('Time(s)');

subplot(3,2,4);
plot(dt*[0:N-1],x1(4,:),'b--');hold on;
plot(dt*[0:N-1],x2(4,:),'r-');
title('$x_4$','Interpreter','latex');
xlabel('Time(s)');

subplot(3,2,5);
plot(dt*[0:N-1],x1(5,:),'b--');hold on;
plot(dt*[0:N-1],x2(5,:),'r-');
title('$x_5$','Interpreter','latex');
xlabel('Time(s)');

subplot(3,2,6);
plot(dt*[0:N-1],x1(6,:),'b--');hold on;
plot(dt*[0:N-1],x2(6,:),'r-');
title('$x_6$','Interpreter','latex');
xlabel('Time(s)');