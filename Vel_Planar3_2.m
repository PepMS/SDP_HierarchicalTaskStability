%% For help type "rtbdemo"
clearvars; clc; close all;

mtTitle = 'IK Stability.: ';

disp(strcat(mtTitle,'Loading libraries'))

addpath(genpath('~/sdpa/share/sdpa/mex'));
addpath(genpath('~/rvctools'));
addpath('functions');
%% Constants definition
d2r = pi/180;

%% Robot creation
L(1) = Link('revolute','d', 0, 'a', 0.5, 'alpha', 0);
L(2) = Link('revolute','d', 0, 'a', 0.3, 'alpha', 0);
L(3) = Link('revolute','d', 0, 'a', 0.2, 'alpha', 0);

robot    = SerialLink(L,      'name', 'Planar_Robot');
%% Task 1 (End-Effector position)
% Desired position
r1d = [0.65;0.17]; % close to a singularity

% Initial position
r10 = [0.7-0.2*cos(0);0.2*sin(0)];

% Gain
L1 = diag([100 100]);
%% Task 2 (End-Effector orientation)
% Desired orientation
angle = -90*d2r;
r2d = angle;

% Initial orientation
angle0 = -70*d2r;
r20 = angle0;

%Gain
L2 = 200;

%% Algorithm

% Initial configuration
R0 = [cos(angle0) -sin(angle0) 0;...
    sin(angle0) cos(angle0) 0;...
    0 0 1];
T0 = [R0 [r10;0]; 0 0 0 1];
q0 = robot.ikunc(T0);
q0 = q0';
q0_d = zeros(robot.n,1);

% Iterators
dt = 0.01;
tt = 0:dt:5;
i = 1;

% Storage variables
QQ   = zeros(robot.n, length(tt));
QQ_d = zeros(robot.n, length(tt));

% Gains
LL1 = zeros(2,length(tt));
LL2 = zeros(1,length(tt));
LL3 = zeros(2,length(tt));

% singular values
L11 = zeros(1,length(tt));
L22 = zeros(1,length(tt));
L12 = zeros(1,length(tt));

% Errors
EE1 = zeros(2,length(tt));
EE2 = zeros(1,length(tt));

PP = zeros(3,3,length(tt));
M  = zeros(3,3,length(tt));

eVAL = zeros(3,length(tt));

% Variables
q = q0;

for t=tt
    
    % Jacobian computation
    J1 = robot.jacob0(q);
    J1(3:6,:)=[];
    
    qSum = sum(q);
    J2 = ones(1,robot.n);
    
    % Gains Calculation
    [L1, L2] = Vel_computeGains_3DOF_2(J1, J2, robot.n);
    % L1 = eye(2)*20;
    % L2 = 20;
    % Null space projectors
    N1 = (eye(robot.n)-pinv(J1)*J1);
    
    % Matrix M construction
    M11 = eye(2);
    M22 = J2*N1*pinv(J2);
    M21 = J2*pinv(J1);
    
    l11 = min(svd(M11));
    l22 = min(svd(M22));
    l21 = max(svd(M21));
       
    % Compute task errors
    r1 = robot.fkine(q);
    r1d_d = [0;0];
    % e1 = r1d - r1.t(1:2);
    e1 = r1d - r1(1:2,4);
    
    r2 = sum(q);
    e2 = r2d - r2;
    
    % Solve CLIK
    % qd = pinv(J1)*(r1d_d + L1*e1) + N1*J2'*L2*e2;
    qd = pinv(J1)*(r1d_d + L1*e1) + pinv(J2*N1)*(L2*e2 - J2*pinv(J1)*L1*e1)
    q = q + qd*dt;
    
    % Store
    QQ(:, i)   = q;   % Joint position
    QQ_d(:, i) = qd;  % Joint velocities
    
    EE1(:, i) = e1; % Error task 1
    EE2(:, i) = e2; % Error task 2
    
    LL1(:,i) = [L1(1) L1(4)]';  % Gain task 1
    LL2(:,i) = L2;              % Gain task 2
    
    L22(i) = min(svd(M22));     % Singular value for task 2 (projected)
    
    M(:,:,i) = [M11*L1, zeros(2,1); ...
        M21*L1, M22*L2];
    
    M_eVAL(:,i) = eig(M(:,:,i));
    
    % Iterators
    i = i + 1;
end

% robot.plot(QQ','fps',50)


%% Plotting tasks

% Plotting task 1
figure;
subplot(2,1,1)
EE_norm = vecnorm(EE1);
plot(tt, EE_norm);
title('Task1 - EE position')
grid on
axis([0 5 -0.1 0.25]);
subplot(2,1,2)
plot(tt, EE2(1,:)*180/pi);
title('Task2 - EE orientation')
axis([0 5 -40 10]);
grid on

%% Plotting joint values
figure
plot(tt, QQ'*180/pi)
title('Joint values')
legend('q1','q2','q3')
grid on
axis([0 5 -150 150])

%% Plotting M eigenvalues
figure
plot(tt, M_eVAL)
axis([0 5 -1 50]);
title('M Eigenvalues')
grid on
%% Plot singular values
figure
plot(tt,L22);
title('Singular values')
grid on
axis([0 5 0 0.15])
%% Tasks gains
figure
subplot(2,1,1)
plot(tt,LL1)
grid on
title('Task Gains');
legend('X pos', 'Y pos');
subplot(2,1,2)
plot(tt,LL2)
grid on