clc
clear all

syms a l m1 m2 I1 I2 g T1 T2 real
syms q1 q2 u1 u2 ud1 ud2 real

%%%%%%% Reference frames %%%%%%%
i = [1 0 0]; j = [0 1 0]; k = [0 0 1];
X1 = sin(q1)*i - cos(q1)*j; Y1 = cos(q1)*i + sin(q1)*j;
X2 = sin(q1+q2)*i - cos(q1+q2)*j; Y2 = cos(q1+q2)*i + sin(q1+q2)*j;

%%%%%%% Position vectors %%%%%%
r_O_G1 = a*X1;
r_O_E = l*X1;
r_E_G2 = a*X2;
r_O_G2 = r_O_E + r_E_G2;
r_O_E2 = r_O_E + l*X2;

%%%%%% Angular velocities and accelerations %%%%%
om1 = u1*k;  om2 = (u1+u2)*k;
al1 = ud1*k; al2 = (ud1+ud2)*k;

%%%%% Velocties (not needed here) and acelerations of masses %%%%
v_O = 0; 
v_G1 = v_O + cross(om1,r_O_G1);
v_E  = v_O + cross(om1,r_O_E);
v_G2 = v_E + cross(om2,r_E_G2);

a_O = 0;
a_G1 = a_O + cross(om1,cross(om1,r_O_G1)) + cross(al1,r_O_G1);
a_E  = a_O + cross(om1,cross(om1,r_O_E))  + cross(al1,r_O_E);
a_G2 = a_E + cross(om2,cross(om2,r_E_G2)) + cross(al2,r_E_G2);

%%%%% Angular Momentum %%%%%%%%%%%%%%%%%%%%%%
M_O =  cross(r_O_G1,-m1*g*j) + cross(r_O_G2,-m2*g*j) + T1*k;
Hdot_O = cross(r_O_G1,m1*a_G1) + I1*al1 + cross(r_O_G2,m2*a_G2) + I2*al2;

M_E = cross(r_E_G2, -m2*g*j) + T2*k;
Hdot_E = cross(r_E_G2, m2*a_G2) + I2*al2;

%%%%% Equations of motion %%%%%%%%%%%%%%%%%%
AMB_O = Hdot_O - M_O;
AMB_E = Hdot_E - M_E;

%%% The k component has the equations of motion %%%
eqn1 = collect(simplify(AMB_O(3)),[ud1 ud2]);
eqn2 = collect(simplify(AMB_E(3)),[ud1 ud2]);
%eqn1 = eqn1 - eqn2; %%There is no need to do this but we do it to make the resulting M matrix symmetric

%%%% We need to write these as [M][alpha] = RHS 
%%%% were M is a 2x2 matrix alpha = [ud1 ud2] 
%%%% and RHS is a 2x1 matrix. 
RHS1 = -subs(eqn1,[ud1 ud2],[0 0]) %negative sign because we want to take RHS to the right eventually
M11 =  subs(eqn1,[ud1 ud2],[1 0]) + RHS1
M12 =  subs(eqn1,[ud1 ud2],[0 1]) + RHS1

RHS2 = -subs(eqn2,[ud1 ud2],[0 0])
M21 =  subs(eqn2,[ud1 ud2],[1 0]) + RHS2
M22 =  subs(eqn2,[ud1 ud2],[0 1]) + RHS2


%%%%%%% Final system [M] [alpha] = [RHS] %%%%%%%%%%%%%%%
M  = [M11 M12; M21 M22];
RHS = [RHS1; RHS2];

%%% Derivation complete %%%
%%% Copy-paste the M and RHS vector to dbpend_drive

%%%% Total energy (as a check on equations) %%%%
KE = simplify(0.5*m1*dot(v_G1,v_G1) + 0.5*m2*dot(v_G2,v_G2)) + 0.5*I1*dot(om1,om1) + 0.5*I2*dot(om2,om2)
PE = m1*g*dot(r_O_G1,j) + m2*g*dot(r_O_G2,j)
r_O_E
r_O_E2

