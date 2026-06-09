function [q1,q2]=inv_kinematics(x,y,a)
  a1=a(1);
  a2=a(2);
  c2=(x*x+y*y-a1*a1-a2*a2)/(2*a1*a2);  %cos(q2)
  s2=sqrt(1-c2*c2); %对应于q2为正号（肘关节朝上，坐标系中逆时针）
%   s2=sqrt(1-c2*c2); %对应于q2为负号（肘关节朝下）
  q2=atan2(s2,c2);
  s1=((a1+a2*cos(q2))*y-a2*sin(q2)*x)/(x*x+y*y);
  c1=((a1+a2*cos(q2))*x+a2*sin(q2)*y)/(x*x+y*y);
  q1=atan2(s1,c1);
end