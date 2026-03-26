clc
clear all

syms a l m1 m2 I1 I2 g T1 T2 I1_1 I1_2 I1_3 I2_1 I2_2 I2_3 I1_12 I1_13 I1_23 I2_12 I2_13 I2_23 real
syms q1 q2 u1 u2 ud1 ud2 real

% I1=[I1_1 0 0;0 I1_2 0;0 0 I1_3];
% I2=[I2_1 0 0;0 I2_2 0;0 0 I2_3];
I1=[I1_1  I1_12 I1_13;
    I1_12 I1_2  I1_23;
    I1_13 I1_23 I1_3];
I2=[I2_1  I2_12 I2_13;
    I2_12 I2_2  I2_23;
    I2_13 I2_23 I2_3];
q = [q1; q2];
u = [u1; u2];
ud = [ud1; ud2];
%%%%%%% Reference frames %%%%%%%
i = [1;0;0]; j = [0;1;0]; k = [0;0;1];
Rzq1 = [cos(q1)*i+sin(q1)*j        -sin(q1)*i+cos(q1)*j                         k];
Rxq2 = [                  i         cos(q2)*j+sin(q2)*k      -sin(q2)*j+cos(q2)*k];
R1 = Rzq1;
R12 = Rzq1*Rxq2;
I1_w = R1*I1*R1';
I2_w = R12*I2*R12';

X1 = R1(:,1);
Y1 = R1(:,2);
Z1 = R1(:,3);
X2 = R12(:,1);
Y2 = R12(:,2);
Z2 = R12(:,3);

% X2_1 = X1
% Y2_1 = cos(q2)*Y1+sin(q2)*Z1
% Z2_1 = -sin(q2)*Y1+cos(q2)*Z1

%%%%%%% Position vectors %%%%%%
r_O_G1 = a*X1;
r_O_E = l*X1;
r_E_G2 = a*Z2;
r_O_G2 = r_O_E + r_E_G2;
r_0_F = r_O_E + l*Z2;

%%%%%% Angular velocities and accelerations %%%%%
om1 = u1*Z1;
om2 = om1 + u2*X2;

%%%%% Velocties (not needed here) and acelerations of masses %%%%
v_O = 0; 
v_G1 = v_O + cross(om1,r_O_G1);
% v_G1_1 = jacobian(r_O_G1,q)*u
% xx=simple(v_G1-v_G1_1)
v_E  = v_O + cross(om1,r_O_E);
v_G2 = v_E + cross(om2,r_E_G2);
% v_G2_1 = jacobian(r_O_G2,q)*u
% xx=simple(v_G2-v_G2_1)

%%%%% Lagrange¡¯s Equations %%%%%%%%%%%%%%%%%%%%%%
KE = 0.5*m1*dot(v_G1,v_G1) + 0.5*m2*dot(v_G2,v_G2) + 0.5*om1'*I1_w*om1 + 0.5*om2'*I2_w*om2;
PE = m1*g*dot(r_O_G1,k) + m2*g*dot(r_O_G2,k);

L = KE - PE;
% L_dq1 = simple(jacobian(L,u1));
% L_dq2 = simple(jacobian(L,u2));
Ld_qd = jacobian(L,u)
% L_dt1 = simple(jacobian(L_dq1,q1)*u1 + jacobian(L_dq1,u1)*ud1+jacobian(L_dq1,q2)*u2 + jacobian(L_dq1,u2)*ud2);
% L_dt2 = simple(jacobian(L_dq2,q1)*u1 + jacobian(L_dq2,u1)*ud1+jacobian(L_dq2,q2)*u2 + jacobian(L_dq2,u2)*ud2);
Ldd_qd_dt = simple(jacobian(Ld_qd,q)*u+jacobian(Ld_qd,u)*ud);
% Ld_q1 = jacobian(L,q1);
% Ld_q2 = jacobian(L,q2);
Ld_q = jacobian(L,q);

L1 = T1-Ldd_qd_dt(1)+Ld_q(1);
L2 = T2-Ldd_qd_dt(2)+Ld_q(2);
%%% The k component has the equations of motion %%%
eqn1 = collect(simple(L1),[ud1 ud2]);
eqn2 = collect(simple(L2),[ud1 ud2]);
%eqn1 = eqn1 - eqn2; %%There is no need to do this but we do it to make the resulting M matrix symmetric

%%%% We need to write these as [M][alpha] = RHS 
%%%% were M is a 2x2 matrix alpha = [ud1 ud2] 
%%%% and RHS is a 2x1 matrix. 
RHS1 = -subs(eqn1,[ud1 ud2],[0 0]); %negative sign because we want to take RHS to the right eventually
M11 =  subs(eqn1,[ud1 ud2],[1 0]) + RHS1;
M12 =  subs(eqn1,[ud1 ud2],[0 1]) + RHS1;

RHS2 = -subs(eqn2,[ud1 ud2],[0 0]);
M21 =  subs(eqn2,[ud1 ud2],[1 0]) + RHS2;
M22 =  subs(eqn2,[ud1 ud2],[0 1]) + RHS2;


%%%%%%% Final system [M] [alpha] = [RHS] %%%%%%%%%%%%%%%
M  = [M11 M12; M21 M22];
RHS = [RHS1; RHS2];

%%% Derivation complete %%%
%%% Copy-paste the M and RHS vector to dbpend_drive

%%%% Total energy (as a check on equations) %%%%
KE = 0.5*m1*dot(v_G1,v_G1) + 0.5*m2*dot(v_G2,v_G2) + 0.5*om1'*I1_w*om1 + 0.5*om2'*I2_w*om2;
PE = m1*g*dot(r_O_G1,k) + m2*g*dot(r_O_G2,k);

r_O_E;
r_0_F;
fid=fopen(   'rhs.m','w');
fprintf(fid, 'function zdot=rhs(t,z,flag,m1, m2, I1, I2, l, a, g)   \n\n');
 
fprintf(fid, 'q1 = z(1);   u1 = z(2);                         \n');
fprintf(fid, 'q2 = z(3);   u2 = z(4);                         \n');
fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3);I2_1 = I2(1,1);I2_2 = I2(2,2);I2_3 = I2(3,3); \n');
fprintf(fid, 'I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);I2_12 = I2(1,2);I2_13 = I2(1,3);I2_23 = I2(2,3); \n');
fprintf(fid, 'T1 = 0; T2 = 0; \n');

fprintf(fid,'M11 = %s; \n', char((M(1,1))) );
fprintf(fid,'M12 = %s; \n', char((M(1,2))) );
fprintf(fid,'\n');
 
fprintf(fid,'M21 = %s; \n', char((M(2,1))) );
fprintf(fid,'M22 = %s; \n', char((M(2,2))) );
fprintf(fid,'\n');

fprintf(fid,'RHS1 = %s; \n', char((RHS(1))) );
fprintf(fid,'RHS2 = %s; \n', char((RHS(2))) );
fprintf(fid,'\n');
 

fprintf(fid,'MM = [M11 M12;                               \n');
fprintf(fid,'     M21 M22];                             \n\n');
 
fprintf(fid,'RHS = [RHS1; RHS2];                      \n\n');
 
fprintf(fid,'X = MM \\ RHS;                                    \n\n');
 
fprintf(fid,'ud1 = X(1);                                       \n');
fprintf(fid,'ud2 = X(2);                                     \n\n');
 
fprintf(fid,'zdot = [u1 ud1 u2 ud2]'';\n'); 
fclose(fid);

fid=fopen(   'energy.m','w');
fprintf(fid, 'function [KE, PE]=energy(t,z,m1, m2, I1, I2, l, a, g)   \n\n');
 
fprintf(fid, 'q1 = z(1);   u1 = z(2);                         \n');
fprintf(fid, 'q2 = z(3);   u2 = z(4);                         \n');
fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3);I2_1 = I2(1,1);I2_2 = I2(2,2);I2_3 = I2(3,3); \n');
fprintf(fid, 'I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);I2_12 = I2(1,2);I2_13 = I2(1,3);I2_23 = I2(2,3); \n');
fprintf(fid,'KE = %s;\n',char(KE)); 
fprintf(fid,'PE = %s;\n',char(PE)); 
fclose(fid);