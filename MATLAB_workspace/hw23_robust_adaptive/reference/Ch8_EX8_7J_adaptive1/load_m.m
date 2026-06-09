%LOAD_M Modify dynamic parameters with load mass.


% mass
  m_pl  = 10; %10 is added  intentionally as load

% inertia moment relative to the center of mass
  I_pl = 0;

% link 2 + load parameters
  % mass [kg]
    m_2 = m_2 + m_pl;

  % inertia first moment
    m2_lC2 = m2_lC2 +11.0; %11 is added  intentionally

  % inertia moment
    I_2 = I_2+ 12.12; %12.12 is added  intentionally

% minimum set of arm + load dynamic parameters
  pi_l(1) = a(1)*m_1 + m1_lC1 + a(1)*m_2;
  pi_l(2) = a(1)*m1_lC1 + I_1 + k_r1^2*I_m1;
  pi_l(3) = a(2)*m_2 + m2_lC2;
  pi_l(4) = a(2)*m2_lC2 + I_2;
  pi_l(5) = I_m2;
