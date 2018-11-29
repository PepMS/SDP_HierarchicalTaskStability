%%
    c = [1 1 1];
    nVars = 3; % Number of varaibles
    nBlock = 5; % Number of LMIs 
    blockStruct = [3 1 1 1 1]; % Positive since it is symmetric
    
    F = cell(nBlock,nVars + 1);
    F{1,1} = zeros(3);
    F{1,2} = [2*l11, -l12, -l13;...
        -l12, 0, 0;...
        -l13, 0 0];
    F{1,3} = [0, 0, 0;...
        0, 2*l22, -l23;...
        0, -l23 0];
    F{1,4} = [0, 0, 0;...
        0, 0, 0;...
        0, 0, 2*l33];
    
    F{2,1} = 25;
    F{2,2} = 1;
    F{2,3} = 0;
    F{2,4} = 0;
    
    F{3,1} = -50;
    F{3,2} = -1;
    F{3,3} = 0;
    F{3,4} = 0;
    
    F{4,1} = -200;
    F{4,2} = 0;
    F{4,3} = -1;
    F{4,4} = 0;
    
    F{5,1} = -200;
    F{5,2} = 0;
    F{5,3} = 0;
    F{5,4} = -1;
    
    OPTION.print='';
    [objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
    
    la = xOpt(1);
    lb = xOpt(2);
    lc = xOpt(3);
    
    %%
    la_ = sdpvar(1,1);
    lb_ = sdpvar(1,1);
    lc_ = sdpvar(1,1);
    
    A1 = [2*l11, -l12, -l13;...
        -l12, 0, 0;...
        -l13, 0 0];
    A2 = [0, 0, 0;...
        0, 2*l22, -l23;...
        0, -l23 0];
    A3 = [0, 0, 0;...
        0, 0, 0;...
        0, 0, 2*l33];
    
    F = [la_==50,A1*la_+A2*lb_+A3*lc_>=0]; 
    %ops = sdpsettings('solver', 'sedumi', 'verbose',0);
    %optimize(F,[],ops);
    optimize(F);
    
    la = double(la_);
    lb = double(lb_);
    lc = double(lc_);
    %%
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
    
    %% Optimizing M backup
        M11 = eye(2);
    M22 = J2*N1*pinv(J2);
    M33 = J3*N12*pinv(J3);
    M21 = J2*pinv(J1);
    M31 = J3*pinv(J1);
    M32 = J3*N1*pinv(J2);
    
    
    
    % L1_1 = sdpvar(1,1);
    L1_ = eye(2)*10;
    L1_2 = sdpvar(1,1);
    L2_  = sdpvar(1,1);
    L3_ = sdpvar(2,2);
    
    M = [M11*L1_, zeros(2,3); ...
        M21*L1_, M22*L2_, zeros(1,2); ...
        M31*L1_, M32*L2_, M33*L3_];
    
    %% YALMIP Sedumi - P
    
    %     L1_ = sdpvar(1,1);
%     L2_ = sdpvar(1,1);
%     L3_ = sdpvar(1,1);
%     lambda = sdpvar(1,1);
%     
%     P = [2*l11*L1_, -l21*L1_, -l31*L1_;...
%         -l21*L1_, 2*l22*L2_, -l32*L2_; ...
%         -l31*L1_, -l32*L2_, 2*l33*L3_];
%       F = [P>=0,L1_==10,lambda*eye(3)-P>=0, L2_>=10];     