global uLINK
lx  = 0.1; ly  = 0.4; lz  = 0.9;         % shape  [m]
mass = 36.0;                             % mass[kg]
MakeRigidBody(1, [lx ly lz], mass);      % Create rectangular data

uLINK(1).p = [0.0, 0.0, 0]'; % Initial position [m]
uLINK(1).R = eye(3);         % Initial posture
uLINK(1).w = [1, 1, 1]';     % I nitial angular velocity [rad/s]
Dtime   = 0.02;           % time step [s]
EndTime = 5.0;            % stop time [s]
time  = 0:Dtime:EndTime;
figure
AX=[-0.5 0.5]; AY=[-0.5 0.5]; AZ=[-0.5 1.0]; % 3D scope range
for n = 1:length(time)
    L = EulerDynamics(1);    % Euler dynamics
    uLINK(1).R = Rodrigues(uLINK(1).w, Dtime) * uLINK(1).R; % Rodrigues
    uLINK(1).w = uLINK(1).w + Dtime * uLINK(1).dw;          % Euler integral
    ShowObject;                 %  display  objcet

end        