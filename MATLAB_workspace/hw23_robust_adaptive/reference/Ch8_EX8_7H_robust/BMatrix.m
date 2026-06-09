%Example for Section 8.7 
%TextBook: Bruno Siciliano, et al. Robotics: Modelling, Planning and Control, Springer,2009

function  B = BMatrix(c_2);
%BMatrix Inertia matrix of two-link planar arm with load mass.
%       B = BMatrix(c_2) returns 2-by-2 inertia matrix of two-link planar arm
%       with load mass, where:
%
%       c_2 = cos(q(2))

global pi_l a k_r2

B(1,1) = a(1)*pi_l(1) + pi_l(2) + (a(2) + 2*a(1)*c_2)*pi_l(3) + pi_l(4);

B(1,2) = (a(2) + a(1)*c_2)*pi_l(3) + pi_l(4) + k_r2*pi_l(5);

B(2,1) = B(1,2);

B(2,2) = a(2)*pi_l(3) + pi_l(4) + k_r2*k_r2*pi_l(5);
