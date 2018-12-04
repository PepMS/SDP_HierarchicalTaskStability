% Developing purposes - Script develop to expand the equations for 
% stability using the Nakamura approach.

clear
syms L1 L2 L3;
syms J1 J2 J3;
syms J1p J2p J3p;
syms J2N1p J3N12p;
syms sig1 sig2 sig3;
syms N1 N12;

q1 = J1p*L1*sig1;
q2 = q1 + J2N1p*(L2*sig2 - J2*q1);
q3 = q2 + J3N12p*(L3*sig3 - J3*q2);

q3_a = J1p*L1*sig1 + ...
    J2N1p*(L2*sig2 - J2*J1p*L1*sig1) + ...
    J3N12p*(L3*sig3 - J3*(J1p*L1*sig1 + J2N1p*(L2*sig2 - J2*J1p*L1*sig1)));

q_sol = expand(q3);
q_expanded = ...
    + J1p*L1*sig1 ...
    - J2*J2N1p*J1p*L1*sig1 ...
    - J3*J3N12p*J1p*L1*sig1 ...
    + J2*J3*J2N1p*J3N12p*J1p*L1*sig1 ...
    + J2N1p*L2*sig2 ...
    - J3*J2N1p*J3N12p*L2*sig2 ...
    + J3N12p*L3*sig3;




