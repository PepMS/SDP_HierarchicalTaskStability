%% For help type "rtbdemo"
clearvars; clc; close all;

mtTitle = 'IK Stability.: ';

disp(strcat(mtTitle,'Loading libraries'))

addpath(genpath('~/sdpa/share/sdpa/mex'));
addpath(genpath('~/rvctools'));
addpath('functions');
addpath('figures');
addpath('plots');
%% Constants definition
d2r = pi/180;

%% Robot creation
L(1) = Link('revolute','d', 0, 'a', 0.5, 'alpha', 0);
L(2) = Link('revolute','d', 0, 'a', 0.3, 'alpha', 0);
L(3) = Link('revolute','d', 0, 'a', 0.2, 'alpha', 0);

robot = SerialLink(L, 'name', 'Planar_Robot');
%% Task 1 (End-Effector position)
% Desired position
r1d = [0.65;0.17]; % close to a singularity
r1d = [0.76; 0.18];

% Initial position
r10 = [0.7-0.2*cos(0);0.2*sin(0)];

% Gain
L1 = diag([100 100]);
%% Task 2 (End-Effector orientation)
% Desired orientation
angle = -90*d2r;
r2d = angle;
r2d = 1.7484;

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
QQct = zeros(robot.n, length(tt));
QQct_d = zeros(robot.n, length(tt));

% Gains
KK1 = zeros(2,length(tt));
KK2 = zeros(1,length(tt));
KK1ct = zeros(2,length(tt));
KK2ct = zeros(1,length(tt));

% singular values
SV1 = zeros(2,length(tt));
SV2 = zeros(1,length(tt));
SV1ct = zeros(2,length(tt));
SV2ct = zeros(1,length(tt));

% Errors
EE1 = zeros(2,length(tt));
EE2 = zeros(1,length(tt));
EE1ct = zeros(2,length(tt));
EE2ct = zeros(1,length(tt));

% M matrix
M   = zeros(3,3,length(tt));
Mct = zeros(3,3,length(tt));

eVAL = zeros(3,length(tt));
eVALct = zeros(3,length(tt));

% Variables
q   = q0;
qct = q0;

for t=tt
    
    % Jacobian computation
    J1   = robot.jacob0(q);
    J1ct = robot.jacob0(qct); 
    J1(3:6,:)   = [];
    J1ct(3:6,:) = [];
    
    qSum   = sum(q);
    qSumct = sum(qct);
    
    J2   = ones(1,robot.n);
    J2ct = ones(1,robot.n);
    
    % Gains Calculation
    [K1, K2] = Vel_computeGains_3DOF_2_Chi(J1, J2, robot.n);
    K1ct = eye(2)*5;
    K2ct = 5;
    
    % Null space projectors
    N1      = (eye(robot.n)-pinv(J1)*J1);
    N1ct    = (eye(robot.n)-pinv(J1ct)*J1ct);
    
    % Matrix M construction
    M11 = eye(2);
    M22 = J2*N1*pinv(J2);
    M21 = J2*pinv(J1);
    M11ct = eye(2);
    M22ct = J2ct*N1ct*pinv(J2ct);
    M21ct = J2ct*pinv(J1ct);
       
    SV1(:, i) = svd(J1);
    SV2(:, i) = svd(N1*pinv(J2));
    SV1ct(:, i) = svd(J1ct);
    SV2ct(:, i) = svd(N1ct*pinv(J2ct));
           
    % Compute task errors
    r1 = robot.fkine(q);
    r1d_d = [0;0];
    %e1 = r1d - r1.t(1:2);
    e1 = r1d - r1(1:2,4);
    r1ct = robot.fkine(qct);
    r1dct_d = [0;0];
    %e1 = r1d - r1.t(1:2);
    e1ct = r1d - r1ct(1:2,4);
   
    r2 = sum(q);
    e2 = r2d - r2;
    r2ct = sum(qct);
    e2ct = r2d - r2ct;
    
    % Solve CLIK
    qd = pinv(J1)*(r1d_d + K1*e1) + N1*pinv(J2)*K2*e2;
    q = q + qd*dt;
    qdct = pinv(J1ct)*(r1dct_d + K1ct*e1ct) + N1ct*pinv(J2ct)*K2ct*e2ct;
    qct = qct + qdct*dt;
    
    % Store
    QQ(:, i)   = q;   % Joint position
    QQ_d(:, i) = qd;  % Joint velocities
    QQct(:, i)   = qct;   % Joint position
    QQct_d(:, i) = qdct;  % Joint velocities
    
    EE1(:, i) = e1; % Error task 1
    EE2(:, i) = e2; % Error task 2
    EE1ct(:, i) = e1ct; % Error task 1
    EE2ct(:, i) = e2ct; % Error task 2
    
    KK1(:,i) = [K1(1) K1(4)]';  % Gain task 1
    KK2(:,i) = K2;              % Gain task 2
    
    M(:,:,i) = [M11*K1, zeros(2,1); ...
        M21*K1, M22*K2];
    
    Mct(:,:,i) = [M11ct*K1, zeros(2,1); ...
        M21*K1, M22*K2];
    
    eVAL(:,i) = eig(M(:,:,i));
    eVALct(:,i) = eig(Mct(:,:,i));
    
    % Iterators
    i = i + 1;
end

%% Plotting robot
conf_num  = size(QQ,2);
conf_num  = 100;
conf_show = 2;
conf_step = cast(conf_num/conf_show,'uint8');

figure(8);
ops =  {'ortho','view','top','noshadow','noshading','notiles','nowrist','noname','jointdiam',5,'linkcolor','g'};
robot.plotopt = ops;
robot.plot(QQ(:,1)');
hold on;
axis([-0.1 1 -0.2 0.6])
j = 1 + conf_step;

%rob_cell = cell(1,conf_show);
for i=1:conf_show
    rob_cell = SerialLink(robot, 'name', strcat('robot',int2str(i)));
    rob_cell.plotopt = ops;
    rob_cell.plot(QQ(:,j)');
    j = j + conf_step;
end
hold off

%% Plotting tasks

% Plotting task 1
task1Err_fig = task1Err_plot(1, tt, [vecnorm(EE1ct);vecnorm(EE1)]);
task2Err_fig = task2Err_plot(2, tt, [EE2ct;EE2]);

genericPrintFig(task1Err_fig, './plots/errorTask1');
genericPrintFig(task2Err_fig, './plots/errorTask2');
%% Plotting joint values
jointVel_fig = jointValues_plot(3, tt, QQ_d);
genericPrintFig(jointVel_fig,'./plots/jointVel');

jointVelct_fig = jointValues_plot(31, tt, QQct_d);

%% Tasks gains
gains_fig = gains_plot(4,tt,[KK1;KK2]);
genericPrintFig(gains_fig,'./plots/gains');

%% Plot singular values
svd_fig = svd_plot(5, tt, [SV1;SV2]);

%% Plotting M eigenvalues

eVal_fig = eVal_plot(6, tt, eVAL);
