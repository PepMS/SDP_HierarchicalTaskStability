function [L1, L2, L3] = computeGains(J1, J2, J3, njoints)

% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

if rank(pinv(J3))+rank(pinv(J12)) ~= rank([pinv(J1) pinv(J2) pinv(J3)])
    disp('Not equaaaal')
end

M11 = eye(2);
M22 = J2*N1*pinv(J2);
M33 = J3*N12*pinv(J3);
M21 = J2*pinv(J1);
M31 = J3*pinv(J1);
M32 = J3*N1*pinv(J2);

M = [M11, zeros(2,3); ...
    M21, M22, zeros(1,2); ...
    M31, M32, M33];

% Find optimal gains - Solving SDP
nVars = 6; % Number of varaibles
nBlock = 8; % Number of LMIs
blockStruct = [5 2 2 1 1 2 2 5]; % Positive since it is symmetric

c = [0 0 0 0 0 1];

% Minimize maximum eigenvalue of M
F = cell(nBlock, nVars + 1);
%F{1,1} = 5*eye(5);
F{1,2} = -M(:,1);
F{1,3} = -M(:,2);
F{1,4} = -M(:,3);
F{1,5} = -M(:,4);
F{1,6} = -M(:,5);
F{1,7} = eye(5);

% Task 1, eigenvalues lower bound
F{2,1} = eye(2)*1;
F{2,2} = [1 0; 0 0];
F{2,3} = [0 0; 0 1];

% Task 1, eigenvalues upperbound
F{3,1} = -eye(2)*100;
F{3,2} = [-1 0; 0 0];
F{3,3} = [0 0; 0 -1];

% Task 2 gains. lower bound
F{4,1} = 1;
F{4,4} = 1;

% Task 2 gains. upper bound
F{5,1} = -100;
F{5,4} = -1;

% Task 3, positive semidefinite (eigen values lower bound)
F{6,1} = eye(2)*1;
F{6,5} = [1 0; 0 0];
F{6,6} = [0 0; 0 1];

% Task 3, eigenvalues upperbound
F{7,1} = -eye(2)*100;
F{7,5} = [-1 0; 0 0];
F{7,6} = [0 0; 0 -1];

% Lower bound eigenvalue of M
%F{8,1} = 0.1*eye(5);
F{8,2} = M(:,1);
F{8,3} = M(:,2);
F{8,4} = M(:,3);
F{8,5} = M(:,4);
F{8,6} = M(:,5);

OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

L1 = double(diag([xOpt(1) xOpt(2)]));
L2 = double(xOpt(3));
L3 = double(diag([xOpt(4) xOpt(5)]));
end