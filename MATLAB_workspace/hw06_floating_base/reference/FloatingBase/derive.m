clc
clear all

syms m g xG yG zG I1_1 I1_2 I1_3 I1_12 I1_13 I1_23 real
syms q1 q2 q3 q4 q5 q6 u1 u2 u3 u4 u5 u6 ud1 ud2 ud3 ud4 ud5 ud6 real
syms lCuboid wCuboid hCuboid real%长方体的长宽高
%浮动基的位置为x=q1,y=q2,z=q3,欧拉姿态角旋转顺序为ZYX，角度为qx=q4, qy=q5,qz=q6
%质量为m，质心在浮动基坐标系中的坐标为xG,yG,zG， 绕质心的惯性张量为I

Ib=[I1_1  I1_12 I1_13;
    I1_12 I1_2  I1_23;
    I1_13 I1_23 I1_3];

q = [q1 q2 q3 q4 q5 q6]';
u = [u1 u2 u3 u4 u5 u6]';
ud = [ud1 ud2 ud3 ud4 ud5 ud6]';
%%%%%%% Reference frames %%%%%%%
i = [1;0;0]; j = [0;1;0]; k = [0;0;1];
Rzqz = Rot('z',q6);
Ryqy = Rot('y',q5);
Rxqx = Rot('x',q4);
R_ZYX = Rzqz*Ryqy*Rxqx;%ZYX欧拉角
Ib_w = R_ZYX*Ib*R_ZYX';
%%%%%%% Position vectors %%%%%%
r_B = [q1 q2 q3]';%基座参考坐标系原点位置
r_B_G = R_ZYX*[xG yG zG]';%
r_G = r_B + r_B_G;%质心位置

p1 = r_G + R_ZYX*[lCuboid wCuboid hCuboid]'/2;%长方体的8个顶点,中心为质心r_G
p2 = r_G + R_ZYX*[-lCuboid wCuboid hCuboid]'/2;
p3 = r_G + R_ZYX*[-lCuboid -wCuboid hCuboid]'/2;
p4 = r_G + R_ZYX*[lCuboid -wCuboid hCuboid]'/2;
p5 = r_G + R_ZYX*[lCuboid wCuboid -hCuboid]'/2;
p6 = r_G + R_ZYX*[-lCuboid wCuboid -hCuboid]'/2;
p7 = r_G + R_ZYX*[-lCuboid -wCuboid -hCuboid]'/2;
p8 = r_G + R_ZYX*[lCuboid -wCuboid -hCuboid]'/2;

r_Bx = r_B + R_ZYX*i;%基座参考坐标系的x、y、z轴单位向量顶点，用于画坐标系
r_By = r_B + R_ZYX*j;
r_Bz = r_B + R_ZYX*k;

r_Gx = r_G + R_ZYX*i;%质心参考坐标系的x、y、z轴单位向量顶点，用于画坐标系
r_Gy = r_G + R_ZYX*j;
r_Gz = r_G + R_ZYX*k;
%%%%%% Angular velocities and accelerations %%%%%
 %om = R_ZYX*(u4*i + Rxqx'*(u5*j)+Rxqx'*Ryqy'*(u6*k));
  om = R_ZYX*u4*i + Rzqz*Ryqy*(u5*j)+Rzqz*(u6*k);
 
 AM = Ib_w*om;
%%%%% Velocties (not needed here) and acelerations of masses %%%%
v_B = [u1 u2 u3]';
v_G = v_B + cross(om,r_B_G);

%%%%% Lagrange’s Equations %%%%%%%%%%%%%%%%%%%%%%
KE = 0.5*m*dot(v_G,v_G) + 0.5*om'*Ib_w*om;
PE = m*g*dot(r_G,k);

L = KE - PE;
% L_dq1 = simple(jacobian(L,u1));
% L_dq2 = simple(jacobian(L,u2));
Ld_qd = jacobian(L,u);
% Ldd_qd_dt = simple(jacobian(Ld_qd,q)*u+jacobian(Ld_qd,u)*ud)
Ldd_qd_dt = jacobian(Ld_qd,q)*u+jacobian(Ld_qd,u)*ud;
% Ld_q1 = jacobian(L,q1);
% Ld_q2 = jacobian(L,q2);
Ld_q = jacobian(L,q);

for ii=1:6
    L(ii) = Ldd_qd_dt(ii) - Ld_q(ii);
    eqn(ii) = collect(L(ii),ud);%%% The k component has the equations of motion %%%
    RHS(ii,1) = -subs(eqn(ii),ud,zeros(6,1)); %negative sign because we want to take RHS to the right eventually
    M(ii,1) =  subs(eqn(ii),ud,[1 0 0 0 0 0]') + RHS(ii,1);%%%% We need to write these as [M][alpha] = RHS 
    M(ii,2) =  subs(eqn(ii),ud,[0 1 0 0 0 0]') + RHS(ii,1);
    M(ii,3) =  subs(eqn(ii),ud,[0 0 1 0 0 0]') + RHS(ii,1);
    M(ii,4) =  subs(eqn(ii),ud,[0 0 0 1 0 0]') + RHS(ii,1);
    M(ii,5) =  subs(eqn(ii),ud,[0 0 0 0 1 0]') + RHS(ii,1);
    M(ii,6) =  subs(eqn(ii),ud,[0 0 0 0 0 1]') + RHS(ii,1);
end
size(RHS)
size(M)

if 1
nDegree=6;

fid=fopen(   'rhs.m','w');

fprintf(fid, 'function zdot=rhs(t,z,flag,m,xG,yG,zG,I1,g)   \n\n');

for ii=1:nDegree
   fprintf(fid, 'q%d = z(%d);   u%d = z(%d);\n',ii,2*ii-1,ii,2*ii);
end
fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3); \n');
fprintf(fid, 'I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);\n\n');
for ii=1:nDegree
    for jj=1:nDegree
        fprintf(fid,'M%d%d =%s;', ii,jj,char(M(ii,jj)));
    end
    fprintf(fid,'\n\n');
end

for ii=1:nDegree 
    fprintf(fid,'RHS%d = %s; \n',ii, char(RHS(ii)) );
end
fprintf(fid,'\n\n');

fprintf(fid,'MM = [');
for ii=1:nDegree
    for jj=1:nDegree
        fprintf(fid,'M%d%d ',ii,jj);
    end
    if ii~=nDegree
        fprintf(fid,';\n');
    end
end
fprintf(fid,' ];\n\n');

fprintf(fid,'RHS = [');
for ii=1:nDegree
    fprintf(fid,'RHS%d',ii);
    if ii~=nDegree
        fprintf(fid,';');
    end
end
fprintf(fid,' ];\n\n');

fprintf(fid,'X = MM \\ RHS;                                    \n\n');

for ii=1:nDegree 
    fprintf(fid,'ud%d = X(%d); \n',ii,ii);
end

fprintf(fid,'zdot = [');
for ii=1:nDegree 
    fprintf(fid,'u%d ud%d ',ii,ii); 
end
fprintf(fid,' ]'';\n\n');
fclose(fid);
% 
fid=fopen(   'energy.m','w');
fprintf(fid, 'function [KE, PE]=energy(t,z,m,xG,yG,zG,I1,g)   \n\n');
 
for ii=1:nDegree
    fprintf(fid, 'q%d = z(%d);   u%d = z(%d);\n',ii,2*ii-1,ii,2*ii);
end
fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3); \n');
fprintf(fid, 'I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);\n\n');

fprintf(fid,'KE = %s;\n',char(KE)); 
fprintf(fid,'PE = %s;\n',char(PE)); 
fclose(fid);

fid=fopen(   'draw_fb.m','w');
fprintf(fid, 'function [r_B, r_Bx, r_By, r_Bz, r_G, r_Gx, r_Gy, r_Gz, p1, p2, p3, p4, p5, p6, p7, p8] = draw_fb(t,z,xG,yG,zG,lCuboid, wCuboid, hCuboid)   \n\n');
for ii=1:nDegree
    fprintf(fid, 'q%d = z(%d);   u%d = z(%d);\n',ii,2*ii-1,ii,2*ii);
end
fprintf(fid,'p1(1,1) = %s;\n',char(p1(1))); 
fprintf(fid,'p1(2,1) = %s;\n',char(p1(2))); 
fprintf(fid,'p1(3,1) = %s;\n',char(p1(3))); 

fprintf(fid,'p2(1,1) = %s;\n',char(p2(1))); 
fprintf(fid,'p2(2,1) = %s;\n',char(p2(2))); 
fprintf(fid,'p2(3,1) = %s;\n',char(p2(3))); 

fprintf(fid,'p3(1,1) = %s;\n',char(p3(1))); 
fprintf(fid,'p3(2,1) = %s;\n',char(p3(2))); 
fprintf(fid,'p3(3,1) = %s;\n',char(p3(3))); 

fprintf(fid,'p4(1,1) = %s;\n',char(p4(1))); 
fprintf(fid,'p4(2,1) = %s;\n',char(p4(2))); 
fprintf(fid,'p4(3,1) = %s;\n',char(p4(3))); 

fprintf(fid,'p5(1,1) = %s;\n',char(p5(1))); 
fprintf(fid,'p5(2,1) = %s;\n',char(p5(2))); 
fprintf(fid,'p5(3,1) = %s;\n',char(p5(3))); 

fprintf(fid,'p6(1,1) = %s;\n',char(p6(1))); 
fprintf(fid,'p6(2,1) = %s;\n',char(p6(2))); 
fprintf(fid,'p6(3,1) = %s;\n',char(p6(3))); 

fprintf(fid,'p7(1,1) = %s;\n',char(p7(1))); 
fprintf(fid,'p7(2,1) = %s;\n',char(p7(2))); 
fprintf(fid,'p7(3,1) = %s;\n',char(p7(3))); 

fprintf(fid,'p8(1,1) = %s;\n',char(p8(1))); 
fprintf(fid,'p8(2,1) = %s;\n',char(p8(2))); 
fprintf(fid,'p8(3,1) = %s;\n',char(p8(3))); 

fprintf(fid,'r_Bx(1,1) = %s;\n',char(r_Bx(1))); 
fprintf(fid,'r_Bx(2,1) = %s;\n',char(r_Bx(2))); 
fprintf(fid,'r_Bx(3,1) = %s;\n',char(r_Bx(3))); 

fprintf(fid,'r_By(1,1) = %s;\n',char(r_By(1))); 
fprintf(fid,'r_By(2,1) = %s;\n',char(r_By(2))); 
fprintf(fid,'r_By(3,1) = %s;\n',char(r_By(3))); 

fprintf(fid,'r_Bz(1,1) = %s;\n',char(r_Bz(1))); 
fprintf(fid,'r_Bz(2,1) = %s;\n',char(r_Bz(2))); 
fprintf(fid,'r_Bz(3,1) = %s;\n',char(r_Bz(3))); 

fprintf(fid,'r_B(1,1) = %s;\n',char(r_B(1))); 
fprintf(fid,'r_B(2,1) = %s;\n',char(r_B(2))); 
fprintf(fid,'r_B(3,1) = %s;\n',char(r_B(3))); 

fprintf(fid,'r_G(1,1) = %s;\n',char(r_G(1))); 
fprintf(fid,'r_G(2,1) = %s;\n',char(r_G(2))); 
fprintf(fid,'r_G(3,1) = %s;\n',char(r_G(3))); 

fprintf(fid,'r_Gx(1,1) = %s;\n',char(r_Gx(1))); 
fprintf(fid,'r_Gx(2,1) = %s;\n',char(r_Gx(2))); 
fprintf(fid,'r_Gx(3,1) = %s;\n',char(r_Gx(3)));

fprintf(fid,'r_Gy(1,1) = %s;\n',char(r_Gy(1))); 
fprintf(fid,'r_Gy(2,1) = %s;\n',char(r_Gy(2))); 
fprintf(fid,'r_Gy(3,1) = %s;\n',char(r_Gy(3)));

fprintf(fid,'r_Gz(1,1) = %s;\n',char(r_Gz(1))); 
fprintf(fid,'r_Gz(2,1) = %s;\n',char(r_Gz(2))); 
fprintf(fid,'r_Gz(3,1) = %s;\n',char(r_Gz(3)));
end

fid=fopen(   'AngularMomentum.m','w');
fprintf(fid, 'function [AMx AMy AMz]=AngularMomentum(t,z,m,xG,yG,zG,I1,g)   \n\n');
 
for ii=1:nDegree
    fprintf(fid, 'q%d = z(%d);   u%d = z(%d);\n',ii,2*ii-1,ii,2*ii);
end
fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3); \n');
fprintf(fid, 'I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);\n\n');

fprintf(fid,'AMx = %s;\n',char(AM(1))); 
fprintf(fid,'AMy = %s;\n',char(AM(2))); 
fprintf(fid,'AMz = %s;\n',char(AM(3))); 
fclose(fid);