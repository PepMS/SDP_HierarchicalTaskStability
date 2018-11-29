% With this function I've fixed the gain related with the first task. 
% When N1*J2+*() degenerates, the gain of the second task is increased 
% in order to keep the error bounded

% Gains of task 1 fixed.
% Optimize for gain 2.
function [L1, L2] = Vel_1_computeGains_3DOF(J1, J2, njoints)

% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

M11 = eye(2)*10;
M22 = J2*N1*pinv(J2);
M21 = J2*pinv(J1)*10*eye(2);

M = [M11, zeros(2,1); ...
    M21, M22];

% Find optimal gains - Solving SDP
nVars = 1; % Number of varaibles
nBlock = 4; % Number of LMIs
blockStruct = [1 1 3 3]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

c = 0;

% Task 2 gains. lower bound
F{1,1} = 1;
F{1,2} = 1;

% Task 2 gains. upper bound
F{2,1} = -1000;
F{2,2} = -1;

% Lower bound eigenvalue of M
F{3,1} = diag([0 0 50])-[M(:,1:2), zeros(3,1)];
F{3,2} = [zeros(3,2), M(:,3)];

%
F{4,1} = diag([0 0 -60])+[M(:,1:2), zeros(3,1)];
F{4,2} = [zeros(3,2), -M(:,3)];

OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

Q = M(:,3)*xOpt(1);
M(:,3) = Q;
eig(M)

L1 = double(diag([10 10]));
L2 = double(xOpt(1));
end