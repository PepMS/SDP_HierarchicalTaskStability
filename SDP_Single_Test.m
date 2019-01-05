% Developing purposes - Script used to test the SDP stability method for
% just one single iteration. 
clear; clc; close all

njoints = 6;

% Load Jacobians
J1 = ...
    [-1.0000   -1.2089   -1.0882   -0.7559   -0.4653   -0.1685 ...
    ; 1.0000    0.5457    0.1330    0.0230    0.0972    0.1414];


J2 = [-0.7660   -0.7660   -0.7660   -0.7660   -0.7660   -0.7660];

J3 = ...
    [-0.2441   -0.4530   -0.3322         0         0         0 ...
    ; 0.9770    0.5228    0.1101         0         0         0];

% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);

M = [J1*pinv(J1) J1*N1*pinv(J2);
    J2*pinv(J1)  J2*N1*pinv(J2)];

la_ = sdpvar(1,1);
lb_ = sdpvar(1,1);
lc_ = sdpvar(1,1);

A1 = [M11,zeros(3,2);M12,0,0;M13,0,0];
A2 = [zeros(3,5);...
zeros(1,3), M22, 0;...
zeros(1,3), M23, 0];
A3 = [zeros(3,5);...
zeros(1,4), M23;...
zeros(1,4), M33];

F = [la_>=25,la_<=50,A1*la_+A2*lb_+A3*lc_>=0]; 
ops = sdpsettings('solver', 'sedumi', 'verbose',0);
%optimize(F,la_+lb_+lc_,ops);
optimize(F,[],ops);

la = double(la_);
lb = double(lb_);
lc = double(lc_);






M11 = eye(2);
M22 = J2*N1*pinv(J2);
M33 = J3*N12*pinv(J3);
M21 = J2*pinv(J1);
M31 = J3*pinv(J1);
M32 = J3*N1*pinv(J2);

L1x = 10;
L1y = 10;
L2_ = 25;
L3x = 150;
L3y = 150;

% M = [M11*diag([L1x L1y]), zeros(2,3); ...
%     M21*diag([L1x L1y]), M22*L2_, zeros(1,2); ...
%     M31*diag([L1x L1y]), M32*L2_, M33*diag([L3x L3y])];

M = [M11, zeros(2,3); ...
    M21, M22, zeros(1,2); ...
    M31, M32, M33];

% Find optimal gains - Solving SDP
nVars = 6; % Number of varaibles
nBlock = 7; % Number of LMIs
blockStruct = [5 2 2 1 1 2 2]; % Positive since it is symmetric

c = [0 0 0 0 0 1];

% Minimize eigenvalues of M
F = cell(nBlock, nVars + 1);
F{1,1} = zeros(0);
F{1,2} = -M(:,1);
F{1,3} = -M(:,2);
F{1,4} = -M(:,3);
F{1,5} = -M(:,4);
F{1,6} = -M(:,5);

% Task 1, positive semidefinite (eigen values lower bound)
F{2,1} = eye(2)*0.1;
F{2,2} = [1 0; 0 0];
F{2,3} = [0 0; 0 1];

% Task 1, eigenvalues upperbound
F{3,1} = -eye(2)*15;
F{3,2} = [-1 0; 0 0];
F{3,3} = [0 0; 0 -1];

% Task 2 gains. lower bound
F{4,1} = 0.1;
F{4,4} = 1;

% Task 2 gains. upper bound
F{5,1} = -30;
F{5,4} = -1;

% Task 3, positive semidefinite (eigen values lower bound)
F{6,1} = eye(2)*0.1;
F{6,5} = [1 0; 0 0];
F{6,6} = [0 0; 0 1];

% Task 3, eigenvalues upperbound
F{7,1} = -eye(2)*200;
F{7,5} = [-1 0; 0 0];
F{7,6} = [0 0; 0 -1];

OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

% Condition 1
disp('Gain task 1, positive definite')
eig(diag([xOpt(1) xOpt(2)]))

L1x = xOpt(1);
L1y = xOpt(2);
L2_ = xOpt(3);
L3x = xOpt(4);
L3y = xOpt(5);

M = [M11*diag([L1x L1y]), zeros(2,3); ...
    M21*diag([L1x L1y]), M22*L2_, zeros(1,2); ...
    M31*diag([L1x L1y]), M32*L2_, M33*diag([L3x L3y])];

% Dot product between matrices (A,B): >> sum(sum(A.*B))
eig(M);