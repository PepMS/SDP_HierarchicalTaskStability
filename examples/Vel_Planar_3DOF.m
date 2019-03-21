%% For help type "rtbdemo"
clearvars; clc; close all;

mtTitle = 'CLIK Stability.: ';

disp(strcat(mtTitle,'Loading libraries'))

addpath(genpath('~/sdpa/share/sdpa/mex'));
addpath(genpath('~/rvctools'));
addpath('../classes');
addpath('../lmi_of');
addpath('../lmi_const');
addpath('../tasks');
addpath('../functions');
addpath('../plots');

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
task_pos = TaskPos_2D(robot, 1, [0.76; 0.18]);
task_pos.bound_u = [300; 300]; % Gain bounds
task_pos.bound_l = [5;5];

task_ori = TaskOri_2D(robot, 2, -70*d2r);
task_ori.bound_u = 2000;
task_ori.bound_l = 1;
%% Task addition to the problem
T = {}; % Cell array of tasks
T{end+1} = task_pos;
T{end+1} = task_ori;

%%% TODO: Method to sort tasks according to their priority
%%% TODO: I am also assuming there is only one task per priority

%% Defining experiment parameters
t_end = 5;
dt = 0.01;

jv_ubound = [10 10 10]';
jv_lbound = [-0.78 -10 -10]';
%% Defining Objective functions and constraints
% Defining which OF we're gonna add
of_LMI = LMI_minMaxEigenvalue();

% LMI_gainBounds = LMI_gainBound(maxGains, minGains);
lmi_gainBounds = LMI_gainBound();
lmi_stability  = LMI_stability();
lmi_jVelBound  = LMI_jVelBound(jv_ubound, jv_lbound);

LMI_l = {};
LMI_l{end+1} = lmi_gainBounds;
LMI_l{end+1} = lmi_stability;
LMI_l{end+1} = lmi_jVelBound;
%% Solve the Gain scheduling problem
clik_SDP = SDPCLIKProblem(robot, q0, T, dt, t_end);

clik_SDP.OF_LMI = of_LMI;
clik_SDP.LMI_l = LMI_l;

clik_SDP = clik_SDP.solve();

%% Plots

t = 0:dt:t_end;
err = vecnorm(clik_SDP.EE).^2;
fig_export = 1;
t1_sat = 0.19;
t2_sat = 0.31;

% Generating figures
g_fig  = plotData(1, 'Gains', t, clik_SDP.KK, '$t$', '$\lambda_', '[s]', '[s$^{-1}$]','');
sv_fig = plotData(2, 'Singular Value', t, clik_SDP.SV(2,:), '$t$', '$\underline{s}_', '[s]', '$\underline{s}_2$', '');
ev_fig = plotData(3, 'Eigenvalues', t, clik_SDP.AA_e, '$t$', '$a_', '[s]', 'EigenValues','');
jv_fig = plotData(4, 'JointVelocities', t, clik_SDP.QQ_d, '$t$', '$\dot{\mbox{\boldmath $q$}}_', '[s]', '[rad/s]','');
% Error plots
ly_fig = plotData(6, 'Lyapunov Function', t, err, '$t$', '', '[s]', '[m$^2$]','log');

% Draw saturation lines
plotHorizontalLine(jv_fig, t, -0.78);
plotVerticalLine(jv_fig, t1_sat);
plotVerticalLine(jv_fig, t2_sat);
plotVerticalLine(ev_fig, t1_sat);
plotVerticalLine(ev_fig, t2_sat);
plotVerticalLine(g_fig, t1_sat);
plotVerticalLine(g_fig, t2_sat);

% Exporting pdf
if fig_export
    genericPrintFig(g_fig,'../plots/gains');
    genericPrintFig(ev_fig,'../plots/AeValues');
    genericPrintFig(jv_fig,'../plots/jointVel');
    genericPrintFig(ly_fig, '../plots/lyapunov');
    genericPrintFig(sv_fig, '../plots/sValues');
end