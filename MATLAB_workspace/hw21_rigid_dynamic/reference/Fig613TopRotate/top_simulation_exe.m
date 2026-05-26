% top_simulation.m
clear
close all
global uLINK G
G = 9.8;
r = 0.2; a = 0.05; c = 0.2;           % 
MakeTop(1, r,a,c);

uLINK(1).p = [0 0 0.3]';                 % initial postion [m]
uLINK(1).R = Rodrigues([1 0 0],pi/50);   % initial posture
uLINK(1).vo= [0 0 0]';                   % initial velocity [m/s]
uLINK(1).w = [0 0 50]';                  % initial angular velocity [rad/s]
Dtime   = 0.002;
EndTime = 2.0;
time  = 0:Dtime:EndTime;
frame_skip = 3;
figure
AX=[-0.2 0.4];  AY=[-0.3 0.3]; AZ=[0 0.8];  %
for n = 1:length(time)
    [f,tau] = TopForce(1);                           % 
    [P,L]   = SE3dynamics(1,f,tau);                  %     
    [uLINK(1).p, uLINK(1).R] = SE3exp(1, Dtime);               
    uLINK(1).w = uLINK(1).w + Dtime * uLINK(1).dw;   
    uLINK(1).vo= uLINK(1).vo+ Dtime * uLINK(1).dvo;  
    if mod(n,frame_skip) == 0
         ShowObject;                                 %
    end
end
