clear; clc; close all

njoints = 6;

% Load Jacobians
J1 = ...
    [-1.0000   -1.2089   -1.0882   -0.7559   -0.4653   -0.1685 ...
    ; 1.0000    0.5457    0.1330    0.0230    0.0972    0.1414];

[U1,E1,V1] = svd(pinv(J1));
[U1,E1,V1] = svd(J1);


J2 = [-0.7660   -0.7660   -0.7660   -0.7660   -0.7660   -0.7660];
[U2,E2,V2] = svd(pinv(J2));
[U2,E2,V2] = svd([zeros(2,6);J2]);

J3 = ...
    [-0.2441   -0.4530   -0.3322         0         0         0 ...
    ; 0.9770    0.5228    0.1101         0         0         0];

sigma1 = [0.6000; -20.2000];
sigma2 = -7.2108;
sigma3 = [6.2106;17.9744];

q1 = pinv(J1)*sigma1;

m1 = 2*(sigma1)'*J1;

nVars = 7; % Number of varaibles
nBlock = 1; % Number of LMIs
blockStruct = 3; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

c = [0 0 0 0 0 0 1];
% Task 2 gains. lower bound
F{1,1} = [zeros(1,blockStruct(1)); zeros(blockStruct(1)-1,1) -eye(blockStruct(1)-1)];
F{1,2} = [m1(1), J1(:,1)'; J1(:,1), zeros(blockStruct(1)-1)];
F{1,3} = [m1(2), J1(:,2)'; J1(:,2), zeros(blockStruct(1)-1)];
F{1,4} = [m1(3), J1(:,3)'; J1(:,3), zeros(blockStruct(1)-1)];
F{1,5} = [m1(4), J1(:,4)'; J1(:,4), zeros(blockStruct(1)-1)];
F{1,6} = [m1(5), J1(:,5)'; J1(:,5), zeros(blockStruct(1)-1)];
F{1,7} = [m1(6), J1(:,6)'; J1(:,6), zeros(blockStruct(1)-1)];
F8 = zeros(blockStruct(1));
F8(1) = 1;
F{1,8} = F8;
% Solve
OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
q1_l = xOpt(1:6);

J2 = eye(6);
nVars = 7; % Number of varaibles
nBlock = 3; % Number of LMIs
blockStruct = [7 2 2]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

c = [0 0 0 0 0 0 1];
% Task 2 gains. lower bound
F{1,1} = [zeros(1,blockStruct(1)); zeros(blockStruct(1)-1,1) -eye(blockStruct(1)-1)];
F{1,2} = [0, J2(:,1)'; J2(:,1), zeros(blockStruct(1)-1)];
F{1,3} = [0, J2(:,2)'; J2(:,2), zeros(blockStruct(1)-1)];
F{1,4} = [0, J2(:,3)'; J2(:,3), zeros(blockStruct(1)-1)];
F{1,5} = [0, J2(:,4)'; J2(:,4), zeros(blockStruct(1)-1)];
F{1,6} = [0, J2(:,5)'; J2(:,5), zeros(blockStruct(1)-1)];
F{1,7} = [0, J2(:,6)'; J2(:,6), zeros(blockStruct(1)-1)];
F8 = zeros(blockStruct(1));
F8(1) = 1;
F{1,8} = F8;

F{2,1} = diag(J1*q1_l);
F{2,2} = diag(J1(:,1));
F{2,3} = diag(J1(:,2));
F{2,4} = diag(J1(:,3));
F{2,5} = diag(J1(:,4));
F{2,6} = diag(J1(:,5));
F{2,7} = diag(J1(:,6));

F{3,1} = -diag(J1*q1_l);
F{3,2} = -diag(J1(:,1));
F{3,3} = -diag(J1(:,2));
F{3,4} = -diag(J1(:,3));
F{3,5} = -diag(J1(:,4));
F{3,6} = -diag(J1(:,5));
F{3,7} = -diag(J1(:,6));

% Solve
OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
q2_l = xOpt(1:6);


e1   = J1*q1   - sigma1
e1_l = J1*q1_l - sigma1





% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);

