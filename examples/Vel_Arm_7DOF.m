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
L(1) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(2) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(3) = Link('revolute','d', 0.3, 'a', 0,    'alpha', -pi/2);
L(4) = Link('revolute','d', 0,   'a', 0.25, 'alpha', -pi/2);
L(5) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(6) = Link('revolute','d', 0,   'a', 0,    'alpha', pi/2);
L(7) = Link('revolute','d', 0.15,'a', 0,    'alpha', 0);

robot = SerialLink(L, 'name', '7DoF_Robot');

%% Robot initial end-effector pose

q0 = d2r*[5 90 0 0 0 90 0]'; %All vectors must be in columns
%% Task Definition
task_pos = TaskPos_3D_EE(robot, 1, [0.40; 0; 0.40]);
task_pos.bound_u = [50; 50; 50]; % Gain bounds
task_pos.bound_l = [50; 50; 50];

task_pos_wr = TaskPos_3D_Wr(robot, 2, 0.10);
task_pos_wr.bound_u = 300;
task_pos_wr.bound_l = 0;

task_pos_el = TaskPos_3D_El(robot, 3, 0.05);
task_pos_el.bound_u = 300;
task_pos_el.bound_l = 0;
%% Task addition to the problem
T = {}; % Cell array of tasks
T{end+1} = task_pos;
T{end+1} = task_pos_wr;
T{end+1} = task_pos_el;

%%% TODO: Method to sort tasks according to their priority
%%% TODO: I am also assuming there is only one task per priority

%% Defining experiment parameters
t_end = 2;
dt = 0.01;

jv_ubound = [ 100  100  100  100  100  100  100]';
jv_lbound = [-100 -100 -100 -100 -100 -100 -100]';
%% Defining Objective functions and constraints
% Defining which OF we're gonna add
of_LMI = LMI_minMaxEigenvalue();

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
norm_err = [vecnorm(clik_SDP.EE(1:3,:)); ...
         abs(clik_SDP.EE(4,:)); ...
         abs(clik_SDP.EE(5,:))];



% Generating figures
g_fig  = plotData(1, 'Gains', t, clik_SDP.KK, '$t$', '$\lambda_', '[s]', '[s$^{-1}$]','');
% Singular Values
ev_fig = plotData(3, 'Eigenvalues', t, clik_SDP.MM_e, '$t$', '$a_', '[s]', 'EigenValues','');
jv_fig = plotData(4, 'JointVelocities', t, clik_SDP.QQ_d, '$t$', '$\dot{\mbox{\boldmath $q$}}_', '[s]', '[rad/s]','');
ee_fig = plotData(5, 'Task Errors', t, norm_err, '$t$', '$e_', '[s]', '[m]','');
ly_fig = plotData(6, 'Lyapunov Function', t, err, '$t$', '', '[s]', '[m$^2$]','log');

% Draw saturation lines
% plotHorizontalLine(jv_fig, t, -0.78);
% plotVerticalLine(jv_fig, t1_sat);
% plotVerticalLine(jv_fig, t2_sat);
% plotVerticalLine(ev_fig, t1_sat);
% plotVerticalLine(ev_fig, t2_sat);

% Exporting pdf
if fig_export
    genericPrintFig(g_fig,'../plots/7dof_gains');
    genericPrintFig(ev_fig,'../plots/7dof_MeValues');
    genericPrintFig(ee_fig,'../plots/7dof_eeValues');
    genericPrintFig(jv_fig,'../plots/7dof_jointVel');
    genericPrintFig(ly_fig, '../plots/7dof_lyapunov');
end