function [L1, L2, L3] = computeGains_bu(J1, J2, J3, njoints)
    
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
    M33 = J3*N12*pinv(J3);
    M21 = J2*pinv(J1);
    M31 = J3*pinv(J1);
    M32 = J3*N1*pinv(J2);
    
    M = [M11, zeros(2,3); ...
        M21, M22, zeros(1,2); ...
        M31, M32, M33];
   
    % Find optimal gains - Solving SDP
    nVars = size(J1,1) + size(J2,1) + size(J3,1); % Number of varaibles
    nBlock = 4; % Number of LMIs 
    blockStruct = [nVars 2 1 2]; % Positive since it is symmetric
        
    c = ones(1,nVars);
    
    F = cell(nBlock, nVars + 1);
    F{1,1} = zeros(nVars);
    F{1,2} = M(:,1);
    F{1,3} = M(:,2);
    F{1,4} = M(:,3);
    F{1,5} = M(:,4);
    F{1,6} = M(:,5);
    
    F{2,1} = 0.1*eye(2);
    F{2,2} = [1 0;0 0];
    F{2,3} = [0 1;0 0];
        
    F{3,1} = 0.1;
    F{3,4} = 1;
    
    F{4,1} = 0.1*eye(2);
    F{4,5} = [1 0;0 0];
    F{4,6} = [0 1;0 0];
        
    OPTION.print='';
    [objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
    
    L1 = diag([xOpt(1) xOpt(2)]);
    L2 = xOpt(3);
    L3 = diag([xOpt(4) xOpt(5)]);

    M_eval = [M11*L1, zeros(2,3); ...
        M21*L1, M22*L2, zeros(1,2); ...
        M31*L1, M32*L2, M33*L3];
    
    eig(M_eval);
end