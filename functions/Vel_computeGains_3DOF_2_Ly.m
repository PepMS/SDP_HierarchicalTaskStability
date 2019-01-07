% With this function all gains are chosen within the optimization framework.
% When N1*J2+*() degenerates, the gain of the second task is increased
% in order to keep the error bounded

function [L1, L2] = Vel_computeGains_3DOF_2_Ly(J1, J2, e1, e2, njoints)

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

A11 = J1*pinv(J1);
A21 = J2*pinv(J1);
A12 = J1*N1*pinv(J2);
A22 = J2*N1*pinv(J2);

A = [A11, A12; ...
    A21, A22];

A_L1 = [A(:,1), zeros(3,2)];
A_L2 = [zeros(3,1), A(:,2), zeros(3,1)];
A_L3 = [zeros(3,2), A(:,3)];

M1 = (A_L1 + A_L1')/2;
M2 = (A_L2 + A_L2')/2;
M3 = (A_L3 + A_L3')/2;

J1p = pinv(J1);
J2p = N1*pinv(J2);
E1 = J1p(:,1)*e1(1);
E2 = J1p(:,2)*e1(2);
E3 = J2p*e2;

% Find optimal gains - Solving SDP
nVars = 4; % Number of varaibles
nBlock = 8; % Number of LMIs
blockStruct = [2 2 1 1 3 3 3 3]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

% Cost function (minimize lambda)
c = [0 0 0 1];

% Task 1 gains. Lower bound
F{1,1} = diag([5 5]);
F{1,2} = [1 0;0 0];
F{1,3} = [0 0;0 1];

% Task 1 gains. Upper bound
F{2,1} = -diag([300 300]);
F{2,2} = [-1 0;0 0];
F{2,3} = [0 0;0 -1];

% Task 2 gains. lower bound
F{3,1} = 1;
F{3,4} = 1;

% Task 2 gains. upper bound
F{4,1} = -2000;
F{4,4} = -1;

% Lower bound eigenvalue of M
F{5,1} = diag([1 1 1]);
F{5,2} = M1;
F{5,3} = M2;
F{5,4} = M3;

% Minimize the maximum eigenvalue of M
F{6,2} = -M1;
F{6,3} = -M2;
F{6,4} = -M3;
F{6,5} = eye(3);

% Upper-Bound for joint velocity
F{7,1} = diag([-10 -10 -10]);
F{7,2} = -diag(E1);
F{7,3} = -diag(E2);
F{7,4} = -diag(E3);

% Lower-Bound for joint velocity
F{8,1} = diag([-0.81 -10 -10]);
F{8,2} = diag(E1);
F{8,3} = diag(E2);
F{8,4} = diag(E3);

% Solve
OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

L1 = double(diag(xOpt(1:2)));
L2 = double(xOpt(3));

A = A*blkdiag(L1,L2);

M = (A+A')/2;


eig(M);
end