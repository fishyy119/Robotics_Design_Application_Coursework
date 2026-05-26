function L = EulerDynamics(j)
%Fig. 6.4 EulerDynamics.m Calculation of Euler¡¯s equation
global uLINK
I = uLINK(j).R * uLINK(j).I * uLINK(j).R';    
L = I * uLINK(j).w;                          
uLINK(j).dw  = I \ (-cross(uLINK(j).w, L));   