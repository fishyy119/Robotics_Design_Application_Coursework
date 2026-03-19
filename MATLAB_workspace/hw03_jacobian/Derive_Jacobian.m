% derivation of jacobian matrix

clear; close all; clc

syms q1 q2 dq1 dq2 'real'  % joint angle and angular rate
syms l1 l2 'real'          % li denotes link length of link i

% unit vectors
i = sym([1; 0]);
j = sym([0; 1]);

% link's unit vector
e1 = sin(q1) * i - cos(q1) * j;
e2 = sin(q1 + q2) * i - cos(q1 + q2) * j;

% joint and end-effector position
P0 = 0 * i + 0 * j;
P1 = P0 + l1 * e1;
P2 = P1 + l2 * e2;

% joint angle
q = [q1; q2];

% matlab jacobian function
J2 = jacobian(P2, q)

dJ_dt = diff(J2, q1)*dq1 + diff(J2, q2)*dq2

%% how derive the jacobian
%  dP2     dP2     dq
% ----- = ----- * -----
%  dt      dq      dt
x = P2(1);
y = P2(2);
dx = diff(x, q1) * dq1 + diff(x, q2) * dq2;
dy = diff(y, q1) * dq1 + diff(y, q2) * dq2;

% velocity of P2
dP2_dt = [dx; dy];

% Jee * dq = dP2_dt
% [Jee, ~] = equationsToMatrix(dP2_dt, [dq1; dq2]);

% test the results
% err = Jee - J2

