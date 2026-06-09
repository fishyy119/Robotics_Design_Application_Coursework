function zdot = DynFunc(t, z, q_ref, dq_ref, ddq_ref, g, pi_m, pi_l, F_v, K_p, K_d, B_h, DQ, rho, ep, a, k_r2)
%DYNFUNC Dynamics function.
%        zdot = DynFunc(t, z, ......) returns 
%        the differential of state vector of system: zdot = [dq1 ddq1 dq2 ddq2]'
%        where: z = [q1 dq1 q2 dq2]' is the state vector and others are
%        parameters which are declared in params.m script.

Q_rel = [z(1), z(3)];
dQ_rel = [z(2), z(4)];

% calculate the errors of joint pos and vel
q_err = q_ref - Q_rel;
dq_err = dq_ref - dQ_rel;
Xi = [q_err dq_err]';

% get robust component w
x = DQ * Xi; % 2*1 matrix
if norm(x) >= ep
    w = rho / norm(x) * x;
else
    w = rho / ep * x;
end

% get the inertia torque
y = w + ddq_ref' + K_p * q_err' + K_d * dq_err';
Inertia_Tau = B_h * y;

% get the friction and gravity compensation
Tau_FGC = F_v * dQ_rel' + [g*pi_m(1)*cos(Q_rel(1)) + g*pi_m(3)*cos(Q_rel(1)+Q_rel(2)); g*pi_m(3)*cos(Q_rel(1)+Q_rel(2))];

% plus Inertia_Tau and Tau_FGC to get total torque
Tau = Inertia_Tau + Tau_FGC;

GG = [g * pi_l(1) * cos(Q_rel(1)) + g * pi_l(3) * cos(Q_rel(1) + Q_rel(2)); g * pi_l(3) * cos(Q_rel(1) + Q_rel(2))];
CC = [ (sin(Q_rel(2)) * (-a(1) * pi_l(3))) * dQ_rel(2) * (2 * dQ_rel(1) + dQ_rel(2)); -(sin(Q_rel(2)) * (-a(1) * pi_l(3))) * dQ_rel(1)^2];

delta_Tau = Tau - F_v * dQ_rel' - GG - CC;

B_q = [a(1)*pi_l(1) + pi_l(2) + (a(2) + 2*a(1)*cos(Q_rel(2)))*pi_l(3) + pi_l(4), (a(2) + a(1)*cos(Q_rel(2)))*pi_l(3) + pi_l(4) + k_r2*pi_l(5);
                   (a(2) + a(1)*cos(Q_rel(2)))*pi_l(3) + pi_l(4) + k_r2*pi_l(5),                  a(2)*pi_l(3) + pi_l(4) + k_r2*k_r2*pi_l(5)];
               
ddq_act = B_q \ delta_Tau;

zdot = [dQ_rel(1) ddq_act(1) dQ_rel(2) ddq_act(2)]';

