function zdot = DynFunc(~, z, q_ref, dq_ref, ddq_ref)
%DYNFUNC Dynamics function.
%        zdot = DynFunc(t, z, q_ref, dq_ref, ddq_ref) returns 
%        the differential of state vector of system: zdot = [dq1 ddq1 dq2 ddq2]'
%        where: z = [q1 dq1 q2 dq2]' is the state vector, 
%        q_ref, dq_ref and ddq_ref are the desired joint position, velocity
%        and acceleration, respectively.

global a k_r2 g F_v pi_l Lam K_ml K_d

q = [z(1) z(3)]';
dq = [z(2) z(4)]';

m_l = z(5);

q_err = q_ref' - q;
dq_err = dq_ref' - dq;

dq_r = Lam * q_err + dq_ref';
ddq_r = Lam * dq_err + ddq_ref';

sig = dq_r - dq;

% obtain joint torque
[Tau, dm_l_nom] = Adaptive_Controller(q, dq, ddq_r, dq_r, sig, m_l);
dm_l = dm_l_nom / K_ml;

% gravitatianal effect
GG = [g * pi_l(1) * cos(q(1)) + g * pi_l(3) * cos(q(1) + q(2)); g * pi_l(3) * cos(q(1) + q(2))];

% centrifugal and Coriolis effect
CC = [ (sin(q(2)) * (-a(1) * pi_l(3))) * dq(2) * (2 * dq(1) + dq(2)); -(sin(q(2)) * (-a(1) * pi_l(3))) * dq(1)^2];

% B(q) * ddq + CC + F_v * dq + GG = Tau + K_d * sig ¡ª¡ª> delta_Tau = B(q) * ddq
delta_Tau = Tau + K_d * sig - F_v * dq - GG - CC;

% B(q) matrix
B_q = [a(1)*pi_l(1) + pi_l(2) + (a(2) + 2*a(1)*cos(q(2)))*pi_l(3) + pi_l(4), (a(2) + a(1)*cos(q(2)))*pi_l(3) + pi_l(4) + k_r2*pi_l(5);
                   (a(2) + a(1)*cos(q(2)))*pi_l(3) + pi_l(4) + k_r2*pi_l(5),                  a(2)*pi_l(3) + pi_l(4) + k_r2*k_r2*pi_l(5)];
               
ddq_act = B_q \ delta_Tau;

zdot = [dq(1) ddq_act(1) dq(2) ddq_act(2) dm_l]';