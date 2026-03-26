function zdot=dbpend_rhs(t,z,flag,m1, m2, I1, ...
                              I2, l, a, g);           

q1 = z(1);                          
u1 = z(2);                          
q2 = z(3);                         
u2 = z(4);                       

w=2;
q1_ref = 1-cos(w*t);
u1_ref = w*sin(w*t);
q2_ref = 1-cos(w*t);
u2_ref = w*sin(w*t);
kp = 1000;
kd=100;
T1 = kp*(q1_ref-q1)+kd*(u1_ref-u1); %zero torques for unactuated system
T2 = kp*(q2_ref-q2)+kd*(u2_ref-u2); 
% T1=0;
% T2=0;
M11 = - I1 - I2 - a^2*m1 - a^2*m2 - l^2*m2 - 2*a*l*m2*cos(q2);%-I1-a^2*m1-m2*l^2;
M12 = - m2*a^2 - l*m2*cos(q2)*a - I2;%-cos(-q1+q2)*l*a*m2;
M21 = - m2*a^2 - l*m2*cos(q2)*a - I2;%-cos(-q1+q2)*l*a*m2;
M22 = - m2*a^2 - I2;%-m2*a^2-I2;

RHS1 = - a*l*m2*sin(q2)*u2^2 - 2*a*l*m2*u1*sin(q2)*u2 - T1 + a*g*m2*sin(q1 + q2) + a*g*m1*sin(q1) + g*l*m2*sin(q1);%m2*g*l*sin(q1)-m2*l*u2^2*a*sin(-q1+q2)+a*sin(q1)*m1*g-T1+T2;
RHS2 = a*l*m2*sin(q2)*u1^2 - T2 + a*g*m2*sin(q1 + q2);%a*sin(q2)*m2*g-T2+m2*a*u1^2*l*sin(-q1+q2);

M    = [M11 M12; M21 M22];
RHS  = [RHS1 ; RHS2];
udot =  M \ RHS;

ud1 = udot(1);
ud2 = udot(2);

zdot = [u1 ud1 u2 ud2]'  ;            
