function [L1, L2, L3] = Vel_computeGains_7DOF_2(J1, J2, J3, njoints)

% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);
% N1 = (eye(njoints)-J1'*J1);
% N12 = (eye(njoints)-J12'*J12);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

if rank(pinv(J3))+rank(pinv(J12)) ~= rank([pinv(J1) pinv(J2) pinv(J3)])
    disp('Not equaaaal')
end

M11 = eye(3);
M22 = J2*N1*pinv(J2);
M33 = J3*N12*pinv(J3);
M21 = J2*pinv(J1);
M31 = J3*pinv(J1);
M32 = J3*N1*pinv(J2);

% M11 = eye(3);
% M22 = J2*N1*J2';
% M33 = J3*N12*J3';
% M21 = J2*J1';
% M31 = J3*J1';
% M32 = J3*N1*J2';

M = [M11, zeros(3,2); ...
    M21, M22, 0; ...
    M31, M32, M33];

% Find optimal gains - Solving SDP
nVars = 6; % Number of varaibles
nBlock = 8; % Number of LMIs
blockStruct = [5 3 3 1 1 1 1 5]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

% Cost function
c = [0 0 0 0 0 1];

% Minimize maximum eigenvalue of M
F{1,2} = [-M(:,1), zeros(5,4)];
F{1,3} = [zeros(5,1), -M(:,2), zeros(5,3)];
F{1,4} = [zeros(5,2), -M(:,3), zeros(5,2)];
F{1,5} = [zeros(5,3), -M(:,4), zeros(5,1)];
F{1,6} = [zeros(5,4), -M(:,5)];
F{1,7} = eye(5);

% Task 1, eigenvalues lowerbound
F{2,1} = eye(3)*50;
F{2,2} = [1 0 0; zeros(2,3)];
F{2,3} = [zeros(1,3); 0 1 0; zeros(1,3)];
F{2,4} = [zeros(2,3); 0 0 1];

% Task 1, eigenvalues upperbound
F{3,1} = -eye(3)*50;
F{3,2} = [-1 0 0; zeros(2,3)];
F{3,3} = [zeros(1,3); 0 -1 0; zeros(1,3)];
F{3,4} = [zeros(2,3); 0 0 -1];

% Task 2 gains. lower bound
F{4,1} = 1;
F{4,5} = 1;

% Task 2 gains. upper bound
F{5,1} = -100;
F{5,5} = -1;

% Task 3, positive semidefinite (eigen values lower bound)
F{6,1} = 1;
F{6,6} = 1;

% Task 3, eigenvalues upperbound
F{7,1} = -100;
F{7,6} = -1;

% Lower bound eigenvalue of M
F{8,1} = diag([0.1 0.1 0.1 0.1 0.1]);
F{8,2} = [M(:,1), zeros(5,4)];
F{8,3} = [zeros(5,1), M(:,2), zeros(5,3)];
F{8,4} = [zeros(5,2), M(:,3), zeros(5,2)];
F{8,5} = [zeros(5,3), M(:,4), zeros(5,1)];
F{8,6} = [zeros(5,4), M(:,5)];

OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

L1 = double(diag([xOpt(1) xOpt(2) xOpt(3)]));
L2 = double(xOpt(4));
L3 = double(xOpt(5));


M11 = M11*L1;
M22 = M22*L2;
M33 = M33*L3;
M21 = M21*L1;
M31 = M31*L1;
M32 = M32*L2;

M = [M11, zeros(3,2); ...
    M21, M22, 0; ...
    M31, M32, M33];

eig(M);
end