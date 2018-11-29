% With this function I've fixed the gain related with the first task. 
% When N1*J2+*() degenerates, the gain of the second task is increased 
% in order to keep the error bounded

% Gains of task 1 fixed.
% Optimize for gain 2.
function [L1, L2] = Vel_computeGains_3DOF_1(J1, J2, njoints)

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

% Set fixed L1 (for this specific solution)
L1 = eye(2)*20;

M11 = eye(2)*L1;
M22 = J2*N1*pinv(J2);
M21 = J2*pinv(J1)*L1;

M = [M11, zeros(2,1); ...
    M21, M22];

% Find optimal gains - Solving SDP
nVars = 2; % Number of varaibles
nBlock = 4; % Number of LMIs
blockStruct = [1 1 3 3]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

c = [0 1];

% Task 2 gains. lower bound
F{1,1} = 1;
F{1,2} = 1;

% Task 2 gains. upper bound
F{2,1} = -1000;
F{2,2} = -1;

% Lower bound eigenvalue of M
F{3,1} = diag([0 0 0.1])-[M(:,1:2), zeros(3,1)];
F{3,2} = [zeros(3,2), M(:,3)];

% Upper bound eigenvalue of M
% F{4,1} = diag([0 0 -100])+[M(:,1:2), zeros(3,1)];
% F{4,2} = [zeros(3,2), -M(:,3)];

% Minimize the maximum eigenvalue of M
F{4,1} = [M(:,1:2), zeros(3,1)];
F{4,2} = [zeros(3,2), -M(:,3)];
F{4,3} = eye(3);

% Solve
OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

L2 = double(xOpt(1));

Q = M(:,3)*L2;
M(:,3)=Q;
eig(M);
end