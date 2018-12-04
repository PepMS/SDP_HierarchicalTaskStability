% With this function all gains are chosen within the optimization framework. 
% When N1*J2+*() degenerates, the gain of the second task is increased 
% in order to keep the error bounded

function [L1, L2] = Vel_computeGains_3DOF_2_Nak(J1, J2, njoints)

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);

if rank(pinv(J1))+rank(pinv(J2)) ~= rank([pinv(J1) pinv(J2)])
    disp('Not equaaaal')
end

M11 = eye(2) - J1*pinv(J2*N1)*J2*pinv(J1);
M12 = J1*pinv(J2*N1);
M22 = J2*pinv(J2*N1);
M21 = J2*(pinv(J1) - pinv(J2*N1)*J2*pinv(J1));

M = [M11, M12; ...
    M21, M22];

% Find optimal gains - Solving SDP
nVars = 4; % Number of varaibles
nBlock = 6; % Number of LMIs
blockStruct = [2 2 1 1 3 3]; % Positive since it is symmetric
F = cell(nBlock, nVars + 1);

% Cost function (minimize lambda)
c = [0 0 0 1];

% Task 1 gains. Lower bound
F{1,1} = diag([5 5]);
F{1,2} = [1 0;0 0];
F{1,3} = [0 0;0 1];

% Task 1 gains. Upper bound
F{2,1} = -diag([1000 1000]);
F{2,2} = [-1 0;0 0];
F{2,3} = [0 0;0 -1];

% Task 2 gains. lower bound
F{3,1} = 1;
F{3,4} = 1;

% Task 2 gains. upper bound
F{4,1} = -1000;
F{4,4} = -1;

% Lower bound eigenvalue of M
F{5,1} = diag([6 6 6]);
F{5,2} = [M(:,1), zeros(3,2)];
F{5,3} = [zeros(3,1), M(:,2), zeros(3,1)];
F{5,4} = [zeros(3,2), M(:,3)];

% Upper bound eigenvalue of M
% F{6,1} = -eye(3)*50;
% F{6,2} = [-M(:,1), zeros(3,2)];
% F{6,3} = [zeros(3,1), -M(:,2), zeros(3,1)];
% F{6,4} = [zeros(3,2), -M(:,3)];


% Minimize the maximum eigenvalue of M
F{6,2} = [-M(:,1), zeros(3,2)];
F{6,3} = [zeros(3,1), -M(:,2), zeros(3,1)];
F{6,4} = [zeros(3,2), -M(:,3)];
F{6,5} = eye(3);

% Solve
OPTION.print = '';
[objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);

L1 = double(diag(xOpt(1:2)));
L2 = double(xOpt(3));

M11 = eye(2)*L1;
M22 = J2*N1*pinv(J2)*L2;
M21 = J2*pinv(J1)*L1;

M = [M11, zeros(2,1); ...
    M21, M22];
eig(M);
end