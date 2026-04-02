clc
clear all

syms a l m1 m2 I1 I2 g T1 T2 I1_1 I1_2 I1_3 I2_1 I2_2 I2_3 real
syms q1 q2 u1 u2 ud1 ud2 real

I1=[I1_1 0 0;0 I1_2 0;0 0 I1_3];
I2=[I2_1 0 0;0 I2_2 0;0 0 I2_3];

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

al1 = jacobian(om1,q)*u + jacobian(om1,u)*ud;
al2 = jacobian(om2,q)*u + jacobian(om2,u)*ud;

%%%%% Velocties (not needed here) and acelerations of masses %%%%
v_O = 0; 
v_G1 = v_O + cross(om1,r_O_G1);
% v_G1_1 = jacobian(r_O_G1,q)*u
% xx=simple(v_G1-v_G1_1)
v_E  = v_O + cross(om1,r_O_E);
v_G2 = v_E + cross(om2,r_E_G2);
% v_G2_1 = jacobian(r_O_G2,q)*u
% xx=simple(v_G2-v_G2_1)

a_O = 0;
a_G1 = a_O + cross(om1,cross(om1,r_O_G1)) + cross(al1,r_O_G1);
% a_G1_1 = jacobian(v_G1,q)*u + jacobian(v_G1,u)*ud
% xx=simple(a_G1-a_G1_1)
a_E  = a_O + cross(om1,cross(om1,r_O_E))  + cross(al1,r_O_E);
a_G2 = a_E + cross(om2,cross(om2,r_E_G2)) + cross(al2,r_E_G2);
% a_G2_1 = jacobian(v_G2,q)*u + jacobian(v_G2,u)*ud
% xx=simple(a_G2-a_G2_1)
%%%%% Angular Momentum %%%%%%%%%%%%%%%%%%%%%%
M_O =  dot(cross(r_O_G1,-m1*g*k) + cross(r_O_G2,-m2*g*k) + T1*Z1,Z1);


Hdot_O = dot(cross(r_O_G1,m1*a_G1) + cross(om1,I1_w*om1) + I1_w*al1 + cross(r_O_G2,m2*a_G2) + cross(om2,I2_w*om2) + I2_w*al2,Z1);

M_E = dot(cross(r_E_G2, -m2*g*k) + T2*X2,X2);
Hdot_E = dot(cross(r_E_G2, m2*a_G2) + cross(om2,I2_w*om2) +I2_w*al2,X2);

%%%%% Equations of motion %%%%%%%%%%%%%%%%%%
AMB_O = M_O - Hdot_O;
AMB_E = M_E - Hdot_E;


%%% The k component has the equations of motion %%%
eqn1 = collect(simple(AMB_O),[ud1 ud2]);
eqn2 = collect(simple(AMB_E),[ud1 ud2]);
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
M  = [M11 M12; M21 M22]
RHS = [RHS1; RHS2]
B = [jacobian(jacobian(RHS1,u1),u2) jacobian(jacobian(RHS2,u1),u2)]'
C = [jacobian(jacobian(RHS1,u1),u1) jacobian(jacobian(RHS1,u2),u2);
     jacobian(jacobian(RHS2,u1),u1) jacobian(jacobian(RHS2,u2),u2);]/2
G = [subs(RHS1,[u1 u2 T1],[0 0 0]) subs(RHS2,[u1 u2 T2],[0 0 0])]'
 

