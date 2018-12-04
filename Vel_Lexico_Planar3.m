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
QQl   = zeros(robot.n, length(tt));
QQl_d = zeros(robot.n, length(tt));

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
EE1l = zeros(2,length(tt));
EE2l = zeros(1,length(tt));

PP = zeros(3,3,length(tt));
M  = zeros(3,3,length(tt));

eVAL = zeros(3,length(tt));

% Variables
q  = q0;
ql = q0;

for t=tt
    
    % Jacobian computation
    J1 = robot.jacob0(q);
    J1(3:6,:)=[];
    
    J1l = robot.jacob0(ql);
    J1l(3:6,:)=[];
    
    qSum = sum(q);
    J2 = ones(1,robot.n);
    
    % Gains Calculation
    L1 = eye(2)*20;
    L2 = 5;
    % Null space projectors
    N1 = (eye(robot.n)-pinv(J1)*J1);
    
    % Matrix M construction
    M11 = eye(2);
    M22 = J2*N1*pinv(J2);
    M21 = J2*pinv(J1);
    
    %L1 = 15*eye(2);
    l11 = min(svd(M11));
    l22 = min(svd(M22));
    l21 = max(svd(M21));
    
    % Compute task errors
    r1 = robot.fkine(q);
    r1d_d = [0;0];
    e1 = r1d - r1(1:2,4);
    
    r2 = sum(q);
    e2 = r2d - r2;
    
    % Solve CLIK
    qd = pinv(J1)*(r1d_d + L1*e1) + N1*J2'*L2*e2;
    q = q + qd*dt;
    
    % ----- Lexicographic ----- %
    r1 = robot.fkine(ql);
    e1l = r1d - r1(1:2,4);
    m1 = 2*(L1*e1l)'*J1l;    
    J1l = [J1l; eye(3)];
    J1l = blkdiag(J1l, -diag(e1l));
    
    r2 = sum(ql);
    e2l = r2d - r2;
    m2 = 2*(L2*e2l)'*J2;
    J2l = [J2; eye(3)];
    
    for kk=1:2
        if kk==1
            nVars = 6; % Number of varaibles
            nBlock = 2; % Number of LMIs
            blockStruct = [8 2]; % Positive since it is symmetric
            F = cell(nBlock, nVars + 1);
            
            % Objective function
            c = [0 0 0 0 0 1];
            
            % Schur Complement
            F1 = [zeros(1,blockStruct(1)); zeros(blockStruct(1)-1,1) -eye(blockStruct(1)-1)];
            F2 = [0, J1l(:,1)'; J1l(:,1), zeros(blockStruct(1)-1)];
            F3 = [0, J1l(:,2)'; J1l(:,2), zeros(blockStruct(1)-1)];
            F4 = [0, J1l(:,3)'; J1l(:,3), zeros(blockStruct(1)-1)];
            F5 = [0, J1l(:,4)'; J1l(:,4), zeros(blockStruct(1)-1)];
            F6 = [0, J1l(:,5)'; J1l(:,5), zeros(blockStruct(1)-1)];
            F7 = zeros(blockStruct(1));
            F7(1) = 1;
            F{1,1} = F1;
            F{1,2} = F2;
            F{1,3} = F3;
            F{1,4} = F4;
            F{1,5} = F5;
            F{1,6} = F6;
            F{1,7} = F7;
            
            % Positive Gains
            F{2,1} = diag([5 5]);
            F{2,5} = [1 0; 0 0];
            F{2,6} = [0 0; 0 1];
            
            % Solve
            OPTION.print = '';
            [objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
            Feval = F2*xOpt(1)+F3*xOpt(2)+F4*xOpt(3)+F5*xOpt(4)+F6*xOpt(5)+F7*xOpt(6) - F1;
            q1_d = xOpt(1:3);
        else
            nVars = 4; % Number of varaibles
            nBlock = 3; % Number of LMIs
            blockStruct = [5 2 2]; % Positive since it is symmetric
            F = cell(nBlock, nVars + 1);
            
            c = [0 0 0 1];
            % Task 2 gains. lower bound
            F{1,1} = [zeros(1,5); zeros(4,1) -eye(4)];
            F{1,2} = [m2(1), J2l(:,1)'; J2l(:,1), zeros(4)];
            F{1,3} = [m2(2), J2l(:,2)'; J2l(:,2), zeros(4)];
            F{1,4} = [m2(3), J2l(:,3)'; J2l(:,3), zeros(4)];
            F5 = zeros(5);
            F5(1) = 1;
            F{1,5} = F5;
            
            F{2,1} = diag(J1l(1:2,:)*q1_d);
            F{2,2} = diag(J1l(1:2,1));
            F{2,3} = diag(J1l(1:2,2));
            F{2,4} = diag(J1l(1:2,3));
            
            F{3,1} = -diag(J1l(1:2,:)*q1_d);
            F{3,2} = -diag(J1l(1:2,1));
            F{3,3} = -diag(J1l(1:2,2));
            F{3,4} = -diag(J1l(1:2,3));
            
            % Solve
            OPTION.print = '';
            [objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
            ql_d = xOpt(1:3);
        end       
    end
    
    % ----- Lexicographic ----- %
        
    % Store
    QQ(:, i)   = q;   % Joint position
    QQ_d(:, i) = qd;  % Joint velocities
    
    ql = ql + ql_d*dt;
    QQl(:, i)   = ql;   % Joint position
    QQl_d(:, i) = ql_d;  % Joint velocities
    
    
    EE1(:, i) = e1; % Error task 1
    EE2(:, i) = e2; % Error task 2
    
    EE1l(:, i) = e1l; % Error task 1
    EE2l(:, i) = e2l; % Error task 2
    
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
% robot.plot(QQl','fps',50)


%% Plotting tasks

% Plotting LS
figure;
subplot(2,1,1)
EE_norm = vecnorm(EE1);
plot(tt, EE_norm);
title('LS: Task1 - EE position')
grid on
axis([0 5 -0.1 0.25]);
subplot(2,1,2)
plot(tt, EE2(1,:)*180/pi);
title('LS: Task2 - EE orientation')
axis([0 5 -40 10]);
grid on

% Plotting Lexico
figure;
subplot(2,1,1)
EE_norm = vecnorm(EE1l);
plot(tt, EE_norm);
title('Lex: Task1 - EE position')
grid on
axis([0 5 -0.1 0.25]);
subplot(2,1,2)
plot(tt, EE2l(1,:)*180/pi);
title('Lex: Task2 - EE orientation')
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
axis([0 5 -1 20]);
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
title('Task Gains');
legend('X pos', 'Y pos');
subplot(2,1,2)
plot(tt,LL2)