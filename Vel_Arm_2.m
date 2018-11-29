%% For help type "rtbdemo"
clearvars; clc; close all;

mtTitle = 'IK Stability.: ';

disp(strcat(mtTitle,'Loading libraries'))

addpath(genpath('~/sdpa/share/sdpa/mex'));
addpath(genpath('~/rvctools'));

%% Constants definition
d2r = pi/180;

%% Robot creation
L(1) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(2) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(3) = Link('revolute','d', 0.3, 'a', 0,    'alpha', -pi/2);
L(4) = Link('revolute','d', 0,   'a', 0.25, 'alpha', -pi/2);
L(5) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(6) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
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

LLa = zeros(3,length(tt));
LLb = zeros(1,length(tt));
LLc = zeros(1,length(tt));

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
M_g  = zeros(5,5,length(tt));

% Variables
q = q0;

for t=tt
    
    % Jacobian computation
    Ja = robot.jacob0(q);
    Ja(4:6,:)=[];
    
    Jb = robot_w.jacob0(q(1:4));
    Jb([1,3:6], :) = [];
    Jb = [Jb, zeros(1,3)];
    
    Jc = robot_e.jacob0(q(1:3));
    Jc([1,3:6], :) = [];
    Jc = [Jc, zeros(1,4)];
    
    [La, Lb, Lc] = Vel_computeGains_7DOF_2(Ja, Jb, Jc, robot.n);
    
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
    
    % Compute task errors
    sa = robot.fkine(q);
    % ea = sad - sa(1:3,end);
    ea = sad - sa.t;
    
    sb = robot_w.fkine(q(1:4));
    % eb = sbd - sb(2,end);
    eb = sbd - sb.t(2);
    
    sc = robot_e.fkine(q(1:3));
    % ec = scd - sc(2,end);
    ec = scd - sc.t(2);
    
    % Solve CLIK
    qd = pinv(Ja)*La*ea + Na*pinv(Jb)*Lb*eb + Nab*pinv(Jc)*Lc*ec;
    q = q + qd*dt;
    
    % Store
    QQ   = [QQ, q];     % Joint position
    QQ_d = [QQ_d, qd];  % Joint velocities
    
    LLa(:, i) = diag(La);   % Gain task a
    LLb(:, i) = Lb;         % Gain task b
    LLc(:, i) = Lc;         % Gain task c
    L11(:, i) = l11;
    L22(:, i) = l22;
    L33(:, i) = l33;
    L12(:, i) = l12;
    L13(:, i) = l13;
    L23(:, i) = l23;
    
    EEa(:, i) = ea;
    EEb(:, i) = norm(eb);
    EEc(:, i) = norm(ec);
    
    M(:,:,i) = [M11, zeros(3,2);...
        M12, M22, M23;...
        M13, M23, M33];
    M_g(:,:,i) = [M11*La, zeros(3,2);...
        M12*La, M22*Lb, M23*Lc;...
        M13*La, M23*Lb, M33*Lc];
    
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
plot(tt,[LLb;LLc])
title('Task b and c gains')
grid on
legend('Task b','Task c')

%% Plot joints

figure
hold on
for ii=1:7
    plot(tt,QQ_d(ii,1:end-1))
end
grid on
title('Joint Velocity')


%% Plot SV
figure
plot(tt,[L22;L33])
title('Minimum singular value of M22 and M33')
legend('L22','L33')
%axis([0 8 0 0.15])
grid on
%% Plotting M eigenvalues
eVAL = zeros(5,length(tt));
for i=1:length(tt)
    eVAL(:,i) = eig(M(:,:,i));
end
figure
plot(tt,eVAL)
title('M eigenvalues');