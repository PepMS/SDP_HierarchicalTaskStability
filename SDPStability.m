%% For help type "rtbdemo"
clearvars; clc; close all;

mtTitle = 'IK Stability.: ';

disp(strcat(mtTitle,'Loading libraries'))

addpath(genpath('~/sdpa/share/sdpa/mex'));
addpath(genpath('~/rvctools'));

%% Constants definition
d2r = pi/180;

%% Robot creation
L(1) = Link('revolute','d', 0,   'a', 0,    'alpha', 90*d2r);
L(2) = Link('revolute','d', 0,   'a', 0,    'alpha', 90*d2r);
L(3) = Link('revolute','d', 0.3, 'a', 0,    'alpha', -90*d2r);
L(4) = Link('revolute','d', 0,   'a', 0.25, 'alpha', -90*d2r);
L(5) = Link('revolute','d', 0,   'a', 0,    'alpha', 90*d2r);
L(6) = Link('revolute','d', 0,   'a', 0,    'alpha', 90*d2r);
L(7) = Link('revolute','d', 0.15,'a', 0,    'alpha', 0);

robot   = SerialLink(L, 'name', 'Planar_Robot');
robot_e = SerialLink(L(1:3),'name','elbow');
robot_w = SerialLink(L(1:4),'name','wrist');
%% Algorithm

% Desired task values
sad = [0.40, 0.0, 0.40]';
sbd = 0.10;
scd = 0.05;

% Initial configuration
q0 = d2r*[5 90 0 0 0 90 0]';
q0_d = zeros(7, 1);

% Iterators
dt = 0.01;
tt = 0:dt:8;
i = 1;

% Storage variables
QQ  = q0;
QQ_d = q0_d;

La = zeros(1,length(tt));
Lb = zeros(1,length(tt));
Lc = zeros(1,length(tt));

L11 = zeros(1,length(tt));
L22 = zeros(1,length(tt));
L33 = zeros(1,length(tt));
L12 = zeros(1,length(tt));
L13 = zeros(1,length(tt));
L23 = zeros(1,length(tt));

EEa = zeros(3,length(tt)); 
EEb = zeros(1,length(tt));
EEc = zeros(1,length(tt));

PP = zeros(3,3,length(tt));
M  = zeros(5,5,length(tt));

% Variables
q = q0;

for t=tt
    
    % Jacobian computation
    Ja = robot.jacob0(q);
    Ja(4:6,:)=[];
        
    Jb = robot_w.jacob0(q(1:4));
    Jb(1, :) = [];
    Jb(2:5, :) = [];
    Jb = [Jb, zeros(1,3)];
        
    Jc = robot_e.jacob0(q(1:3));
    Jc(1, :) = [];
    Jc(2:5, :) = [];
    Jc = [Jc, zeros(1,4)];
    
    % Augmented Jacobians
    Jab = [Ja;Jb];
    
    % Null space projectors
    Na = (eye(7)-pinv(Ja)*Ja);
    Nab = (eye(7)-pinv(Jab)*Jab);
    
    if rank(pinv(Ja))+rank(pinv(Jb)) ~= rank([pinv(Ja) pinv(Jb)])
        disp('Not equaaaal')
    elseif rank(Jc')+rank(Jab') ~= rank([Jc' Jab'])
        disp('Not equaaaal')
    end
        
    % Singular Values
    M11 = eye(3);
    M22 = Jb*Na*pinv(Jb);
    M33 = Jc*Nab*pinv(Jc);
    M12 = Jb*pinv(Ja);
    M13 = Jc*pinv(Ja);
    M23 = Jc*Na*pinv(Jb);
    
    l11 = min(svd(M11));
    l22 = min(svd(M22));
    l33 = min(svd(M33));
    l12 = max(svd(M12));
    l13 = max(svd(M13));
    l23 = max(svd(M23));
    
    % Find optimal gains - Solving SDP
    c = [1 1 1];
    nVars = 3; % Number of varaibles
    nBlock = 5; % Number of LMIs 
    blockStruct = [3 1 1 1 1]; % Positive since it is symmetric
    
    F = cell(nBlock,nVars + 1);
    F{1,1} = 3*eye(3);
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
    
    F{4,1} = -5000;
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
    
    % Compute task errors
    sa = robot.fkine(q);
    ea = sad - sa.t;
    
    sb = robot_w.fkine(q(1:4));
    eb = sbd - sb.t(2);
    
    sc = robot_e.fkine(q(1:3));
    ec = scd - sc.t(2);
    
    % Solve CLIK
    qd = pinv(Ja)*eye(3)*la*ea + Na*pinv(Jb)*lb*eb + Nab*pinv(Jc)*lc*ec;
    q = q + qd*dt;
    
    % Store
    QQ   = [QQ, q];     % Joint position
    QQ_d = [QQ_d, qd];  % Joint velocities
    
    La(:, i) = la;      % Gain task a    
    Lb(:, i) = lb;      % Gain task b
    Lc(:, i) = lc;      % Gain task c
    L11(:, i) = l11;
    L22(:, i) = l22;
    L33(:, i) = l33;
    L12(:, i) = l12;
    L13(:, i) = l13;
    L23(:, i) = l23;
        
    EEa(:, i) = ea;
    EEb(:, i) = norm(eb);
    EEc(:, i) = norm(ec);
    
    PP(:,:,i) = [2*l11*la, -l12*la, -l13*la;...
                -l12*la,    2*l22*lb, -l23*lb;...
                -l13*la, -l23*lb, 2*l33*lc];
    
            %eig(PP(:,:,i))
    M(:,:,i) = [M11*la, zeros(3,2);...
            M12*la, M22*lb, M23*lc;...
            M13*la, M23*lb, M33*lc];
    
    % Iterators
    i = i + 1;
    disp(i)
end

% robot.plot(QQ','fps',250)

%% Plotting Lyapunov Function

% Plotting task 1
figure;
error = [EEa;EEb;EEc];
semilogy(tt, vecnorm(error).^2);
title('Lyapunov function')
grid on


%% Plotting tasks

figure;
plot(tt, vecnorm(EEa),tt,abs(EEb),tt,abs(EEc));
grid on
legend('a','b','c');
title('Tasks errors')

%% Plotting 
figure
plot(tt,La,tt,Lb,tt,Lc)
title('Lambdas ii')
grid on

%% Plotting gains

figure
plot(tt, Lb,tt,Lc)
axis([0 8 0 300])
title('Lambda b and c')
grid on

%% Plot joints

figure
hold on
for ii=1:7
    plot(tt,QQ_d(ii,1:end-1))
end
grid on
title('Joint Velocity')
