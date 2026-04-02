function [KE, PE] = dbpend_energy(t,z,m1, m2, I1, I2, l, a, g)

q1 = z(1);                          
u1 = z(2);                          
q2 = z(3);                         
u2 = z(4);  

KE = (m2*((l*u1*cos(q1) + a*cos(q1 + q2)*(u1 + u2))^2 + (l*u1*sin(q1) + a*sin(q1 + q2)*(u1 + u2))^2))/2 + (I2*(u1 + u2)^2)/2 + (I1*u1^2)/2 + (a^2*m1*u1^2)/2;%(m1*a^2*u1^2)/2 + (m2*a^2*u2^2)/2 + m2*cos(q2)*a*l*u1*u2 + (m2*l^2*u1^2)/2 + (I1*u1^2)/2 + (I2*u2^2)/2;
 
PE = - g*m2*(a*cos(q1 + q2) + l*cos(q1)) - a*g*m1*cos(q1);%-a*m1*g*cos(q1)-m2*g*(l*cos(q1)+a*cos(q2));