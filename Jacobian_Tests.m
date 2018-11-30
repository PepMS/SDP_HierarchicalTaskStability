clear; clc; close all

njoints = 6;

% Load Jacobians
J1 = ...
    [-1.0000   -1.2089   -1.0882   -0.7559   -0.4653   -0.1685 ...
    ; 1.0000    0.5457    0.1330    0.0230    0.0972    0.1414];
[U1,E1,V1] = svd(pinv(J1));


J2 = [-0.7660   -0.7660   -0.7660   -0.7660   -0.7660   -0.7660];
[U2,E2,V2] = svd(pinv(J2));

J3 = ...
    [-0.2441   -0.4530   -0.3322         0         0         0 ...
    ; 0.9770    0.5228    0.1101         0         0         0];

sigma1 = [0.6000; -20.2000];
sigma2 = -7.2108;
sigma3 = [6.2106;17.9744];

q1 = pinv(J1)*sigma1;
e1 = J1*q1 - sigma1;

q2 = pinv(J2)*sigma2;
e2 = J2*q2 - sigma2;
q2_proj = N1*q2;
e2_real = J2*q2_proj-sigma2;
q12 = q1+q2_proj;

J1*q12-sigma1
J2*q12-sigma2

% Augmented Jacobians
J12 = [J1;J2];

% Null space projectors
N1 = (eye(njoints)-pinv(J1)*J1);
N12 = (eye(njoints)-pinv(J12)*J12);

