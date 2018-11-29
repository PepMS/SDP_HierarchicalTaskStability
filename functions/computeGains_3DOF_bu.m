function [L1, L2] = computeGains_3DOF(J1, J2, njoints)

% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

M11 = eye(2);
M22 = J2*N1*pinv(J2);
M21 = J2*pinv(J1);

M = [M11, zeros(2,1); ...
    M21, M22];

% Find optimal gains - Solving SDP
nVars = 4; % Number of varaibles
nBlock = 6; % Number of LMIs
blockStruct = [3 2 2 1 1 3]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);


c = [0 0 0 0];

% Minimize maximum eigenvalue of M
F{6,1} = -100*eye(3);
F{6,2} = -M(:,1);
F{6,3} = -M(:,2);
F{6,4} = -M(:,3);

% Task 1, eigenvalues lower bound
F{2,1} = eye(2)*15;
F{2,2} = [1 0; 0 0];
F{2,3} = [0 0; 0 1];

% Task 1, eigenvalues upperbound
F{3,1} = -eye(2)*15;
F{3,2} = [-1 0; 0 0];
F{3,3} = [0 0; 0 -1];

% Task 2 gains. lower bound
F{4,1} = 1;
F{4,4} = 1;

% Task 2 gains. upper bound
F{5,1} = -1000;
F{5,4} = -1;

% Lower bound eigenvalue of M
F{6,1} = 20*eye(3);
F{6,2} = M(:,1);
F{6,3} = M(:,2);
F{6,4} = M(:,3);

OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

L1 = double(diag([xOpt(1) xOpt(2)]));
L2 = double(xOpt(3));
end