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
L(1) = Link('revolute','d', 0, 'a', 0.50, 'alpha', 0);
L(2) = Link('revolute','d', 0, 'a', 0.43, 'alpha', 0);
L(3) = Link('revolute','d', 0, 'a', 0.35, 'alpha', 0);
L(4) = Link('revolute','d', 0, 'a', 0.30, 'alpha', 0);
L(5) = Link('revolute','d', 0, 'a', 0.30, 'alpha', 0);
L(6) = Link('revolute','d', 0, 'a', 0.22, 'alpha', 0);

robot    = SerialLink(L,      'name', 'Planar_Robot');
robot_3L = SerialLink(L(1:3), 'name', '3rd Link'); 
%% Task 1 (End-Effector position)
% Desired position
r1d = [1.6;0.6];

% Initial position
r10 = [1;1];

% Gain
L1 = diag([10 10]);
%% Task 2 (End-Effector orientation)
% Desired orientation
angle = 60;
r2d = cosd(angle);

% Initial orientation
angle0 = 50;
r20 = cosd(angle0);

%Gain
L2 = 25;

%% Task 3 (Link 3 position)
% Desired position
r3d = [1.1;0.6];

% Initial position
r30 = [0.7;0.6];

% Gain
L3 = diag([150 150]);

%% Algorithm

% Initial configuration
R0 = [cosd(angle0) -sind(angle0) 0;...
    sind(angle0) cosd(angle0) 0;...
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
L33 = zeros(2,length(tt));
L12 = zeros(1,length(tt));
L13 = zeros(1,length(tt));
L23 = zeros(1,length(tt));

% Errors
EE1 = zeros(2,length(tt)); 
EE2 = zeros(1,length(tt));
EE3 = zeros(2,length(tt));

PP = zeros(3,3,length(tt));
M  = zeros(5,5,length(tt));

eVAL = zeros(5,length(tt));

% Variables
q = q0;

for t=tt
    
    % Jacobian computation
    J1 = robot.jacob0(q);
    J1(3:6,:)=[];
        
    qSum = sum(q);
    J2 = -sin(qSum)*ones(1,robot.n);
        
    J3 = robot_3L.jacob0(q(1:3));
    J3(3:6, :) = [];
    J3 = [J3, zeros(2,robot.n-3)];
    
    % Gains Calculation
    [L1, L2, L3] = computeGains(J1, J2, J3, robot.n);
    
    % Augmented Jacobians
    J12 = [J1;J2];
    
    % Null space projectors
    N1 = (eye(robot.n)-pinv(J1)*J1);
    N12 = (eye(robot.n)-pinv(J12)*J12);
    
    % Matrix M construction
    M11 = eye(2);
    M22 = J2*N1*pinv(J2);
    M33 = J3*N12*pinv(J3);
    M21 = J2*pinv(J1);
    M31 = J3*pinv(J1);
    M32 = J3*N1*pinv(J2);
    
    % Compute task errors
    r1 = robot.fkine(q);
    e1 = r1d - r1(1:2,4);
    
    r2 = cos(sum(q));
    e2 = r2d - r2;
    
    r3 = robot_3L.fkine(q(1:3));
    e3 = r3d - r3(1:2,4);
    
    % Solve CLIK
    qd = pinv(J1)*L1*e1 + N1*pinv(J2)*L2*e2 + N12*pinv(J3)*L3*e3;
    % qd = pinv([J1; J2; J3])*[e1;e2;e3];
    q = q + qd*dt;
    
    % Store
    QQ(:, i)   = q;     % Joint position
    QQ_d(:, i) = qd;  % Joint velocities
           
    EE1(:, i) = e1;
    EE2(:, i) = e2;
    EE3(:, i) = e3;
    
    LL1(:,i) = [L1(1) L1(4)]';
    LL2(:,i) = L2;
    LL3(:,i) = [L3(1) L3(4)]';
    
    L22(i) = min(svd(M22));
    L33(:,i) = svd(M33);
    
    M(:,:,i) = [M11*L1, zeros(2,3); ...
        M21*L1, M22*L2, zeros(1,2); ...
        M31*L1, M32*L2, M33*L3];
    eVAL(:,i) = eig(M(:,:,i));
    
    % Iterators
    i = i + 1;
end

% robot.plot(QQ','fps',50)


%% Plotting tasks

% Plotting task 1
figure;
subplot(3,1,1)
plot(tt, EE1(1,:)); hold on
plot(tt, EE1(2,:));
title('Task1 - EE position')
grid on
legend('Xpos','Ypos')
subplot(3,1,2)
plot(tt, EE2(1,:));
title('Task2 - EE orientation')
grid on

subplot(3,1,3)
plot(tt, EE3(1,:)); hold on
plot(tt, EE3(2,:))
title('Task3 - Link 3 position')
grid on
legend('Xpos','Ypos')
%% Plotting joint values
figure
plot(tt, QQ')
title('Joint values')
legend('q1','q2','q3','q4','q5','q6')
grid on

%% Plotting M eigenvalues
figure
plot(tt, eVAL)
axis([0 5 -1 150]);
title('M Eigenvalues')
grid on

%% Plot singular values

figure
plot(tt,[L22;L33(2,:)./L33(1,:)]');
title('Singular values')
grid on
legend('lambda22', 'lambda33');
axis([0 5 -1 0.5])
%% Tasks gains
figure
subplot(3,1,1)
plot(tt,LL1)
title('Task Gains');
subplot(3,1,2)
plot(tt,LL2)
subplot(3,1,3)
plot(tt,LL3)
%% Robot plot
% figure
% hold on
% for i=1:length(tt)
%     robot.plot(QQ(:,i)');
%     text = strcat('iteration',i);
%     title(text)
%     w = waitforbuttonpress
% end

