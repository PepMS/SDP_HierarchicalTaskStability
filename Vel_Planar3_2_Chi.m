%% For help type "rtbdemo"
clearvars; clc; close all;

mtTitle = 'IK Stability.: ';

disp(strcat(mtTitle,'Loading libraries'))

addpath(genpath('~/sdpa/share/sdpa/mex'));
addpath(genpath('~/rvctools'));
addpath('functions');
addpath('classes');
addpath('figures');
addpath('plots');
%% Constants definition
d2r = pi/180;

%% Robot creation
L(1) = Link('revolute','d', 0, 'a', 0.5, 'alpha', 0);
L(2) = Link('revolute','d', 0, 'a', 0.3, 'alpha', 0);
L(3) = Link('revolute','d', 0, 'a', 0.2, 'alpha', 0);

robot = SerialLink(L, 'name', 'Planar_Robot');

%% Robot initial end-effector pose
pos0 = [0.7-0.2*cos(0);0.2*sin(0)];
ori0 = -65*d2r;

q0 = robot.ikunc(from2DPose2T(pos0,ori0));
q0 = q0'; %All vectors must be in columns

%% Task Definition
task_pos = TaskPos(robot, 1, [0.76; 0.18]);
task_ori = TaskOri(robot, 2, -70*d2r);

%% Task addition to the problem
T = {}; % Cell array of tasks
T{end+1} = task_pos;
T{end+1} = task_ori;

%%% TODO: Method to sort tasks according to their priority

%% Algorithm
q0_d = zeros(robot.n, 1);

% Iterators
dt = 0.01;
tt = 0:dt:5;
i = 1;

% Storage variables
QQ   = zeros(robot.n, length(tt));
QQ_d = zeros(robot.n, length(tt));
QQct = zeros(robot.n, length(tt));
QQct_d = zeros(robot.n, length(tt));

% Gains
KK1 = zeros(2,length(tt));
KK2 = zeros(1,length(tt));
KK1ct = zeros(2,length(tt));
KK2ct = zeros(1,length(tt));

% singular values
SV1 = zeros(1,length(tt));
SV2 = zeros(1,length(tt));
SV1ct = zeros(1,length(tt));
SV2ct = zeros(1,length(tt));

% Errors
EE1 = zeros(2,length(tt));
EE2 = zeros(1,length(tt));
EE1ct = zeros(2,length(tt));
EE2ct = zeros(1,length(tt));

% M matrix
M   = zeros(3,3,length(tt));
Mct = zeros(3,3,length(tt));

MeVAL = zeros(3,length(tt));
AeVAL = zeros(3,length(tt));

McteVAL = zeros(3,length(tt));

% Variables
q   = q0;
qct = q0;

% Define a problem variable
clik_SDP = SDP_CLIK_GainScheduling(robot, T, tt);
clik_SDP.solve();



%% Plotting robot
conf_num  = size(QQ,2);
conf_num  = 100;
conf_show = 2;
conf_step = cast(conf_num/conf_show,'uint8');

figure(8);
ops =  {'ortho','view','top','noshadow','noshading','notiles','nowrist','noname',...
    'jointcolor','k','jointdiam',2,'linkcolor','k'};
robot.plotopt = ops;
robot.plot(QQ(:,1)');
%robot.plot(QQ');
%hold on;
axis([-0.1 1 -0.2 0.6])
j = 1 + conf_step;

%rob_cell = cell(1,conf_show);
% for i=1:conf_show
%     rob_cell = SerialLink(robot, 'name', strcat('robot',int2str(i)));
%     rob_cell.plotopt = ops;
%     rob_cell.plot(QQ(:,j)');
%     j = j + conf_step;
% end
%hold off

%% Plotting tasks

% Plotting task 1
task1Err_fig = task1Err_plot(1, tt, [vecnorm(EE1ct);vecnorm(EE1)]);
task2Err_fig = task2Err_plot(2, tt, [EE2ct;EE2]);

genericPrintFig(task1Err_fig, './plots/errorTask1');
genericPrintFig(task2Err_fig, './plots/errorTask2');
%% Plotting joint values
jointVel_fig = jointValues_plot(3, tt, QQ_d);
genericPrintFig(jointVel_fig,'./plots/jointVel');

% jointVelct_fig = jointValues_plot(31, tt, QQct_d);

%% Tasks gains
gains_fig = gains_plot(4,tt,[KK1;KK2]);
genericPrintFig(gains_fig,'./plots/gains');

%% Plot singular values
svd_fig = svd_plot(5, tt, [SV1;SV2]);
genericPrintFig(svd_fig,'./plots/sValues');
%% Plotting A eigenvalues
MeVal_fig = eVal_plot(6, tt, MeVAL);
genericPrintFig(MeVal_fig,'./plots/MeValues');

%% Plotting A eigenvalues
AeVal_fig = eVal_plot(7, tt, AeVAL);
genericPrintFig(AeVal_fig,'./plots/AeValues');

%% Plot Lyapunov function

%% Plotting Lyapunov Function

% Plotting task 1
% figure;
% error = [EEa;EEb;EEc];
% semilogy(tt, vecnorm(error).^2);
% title('Lyapunov function')
% grid on