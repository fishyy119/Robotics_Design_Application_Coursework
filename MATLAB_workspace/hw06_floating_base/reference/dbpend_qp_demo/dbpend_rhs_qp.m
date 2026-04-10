   function zdot=dbpend_rhs_qp(t,z,flag,m1, m2, I1, ...
                              I2, l, a, g, T)           

q1 = z(1);                          
u1 = z(2);                          
q2 = z(3);                        
u2 = z(4);                       

% [M][alpha] = RHS + T
M11 = I1 + I2 + a^2*m1 + a^2*m2 + l^2*m2 + 2*a*l*m2*cos(q2);%-I1-a^2*m1-m2*l^2;
M12 = m2*a^2 + l*m2*cos(q2)*a + I2;%-cos(-q1+q2)*l*a*m2;
M21 = m2*a^2 + l*m2*cos(q2)*a + I2;%-cos(-q1+q2)*l*a*m2;
M22 = m2*a^2 + I2;%-m2*a^2-I2;
RHS1 = a*l*m2*sin(q2)*u2^2 + 2*a*l*m2*u1*sin(q2)*u2 - a*g*m2*sin(q1 + q2) - a*g*m1*sin(q1) - g*l*m2*sin(q1);%m2*g*l*sin(q1)-m2*l*u2^2*a*sin(-q1+q2)+a*sin(q1)*m1*g-T1+T2;
RHS2 = - a*l*m2*sin(q2)*u1^2 - a*g*m2*sin(q1 + q2);%a*sin(q2)*m2*g-T2+m2*a*u1^2*l*sin(-q1+q2);
M    = [M11 M12; M21 M22];
RHS  = [RHS1;RHS2];

% alpha = [M^-1][RHS + T]
udot =  M \ [RHS+T];
ud1 = udot(1);
ud2 = udot(2);

zdot = [u1 ud1 u2 ud2]'  ;     
