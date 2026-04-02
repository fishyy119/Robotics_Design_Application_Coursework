function zdot=rhs(t,z,flag,m1, m2, I1, I2, l, a, g)   

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3);I2_1 = I2(1,1);I2_2 = I2(2,2);I2_3 = I2(3,3); 
I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);I2_12 = I2(1,2);I2_13 = I2(1,3);I2_23 = I2(2,3); 
T1 = 0; T2 = 0; 
M11 = (I2_2*cos(2*q2))/2 - I2_2/2 - I2_3/2 - a^2*m1 - (a^2*m2)/2 - l^2*m2 - I1_3 - (I2_3*cos(2*q2))/2 - I2_23*sin(2*q2) + (a^2*m2*cos(2*q2))/2; 
M12 = a*l*m2*cos(q2) - I2_12*sin(q2) - I2_13*cos(q2); 

M21 = a*l*m2*cos(q2) - I2_12*sin(q2) - I2_13*cos(q2); 
M22 = - I2_1 - a^2*m2; 

RHS1 = I2_12*u2^2*cos(q2) - T1 - I2_13*u2^2*sin(q2) + 2*I2_23*u1*u2*cos(2*q2) + I2_2*u1*u2*sin(2*q2) - I2_3*u1*u2*sin(2*q2) + a*l*m2*u2^2*sin(q2) + a^2*m2*u1*u2*sin(2*q2); 
RHS2 = (I2_3*u1^2*sin(2*q2))/2 - I2_23*u1^2*cos(2*q2) - (I2_2*u1^2*sin(2*q2))/2 - T2 - a*g*m2*sin(q2) - (a^2*m2*u1^2*sin(2*q2))/2; 

MM = [M11 M12;                               
     M21 M22];                             

RHS = [RHS1; RHS2];                      

X = MM \ RHS;                                    

ud1 = X(1);                                       
ud2 = X(2);                                     

zdot = [u1 ud1 u2 ud2]';
