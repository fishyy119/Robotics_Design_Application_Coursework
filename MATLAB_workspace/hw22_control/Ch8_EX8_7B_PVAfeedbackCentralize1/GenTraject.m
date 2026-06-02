%Example Trajectory for Section 8.7 
%TextBook: Bruno Siciliano, et al. Robotics: Modelling, Planning and Control, Springer,2009


%--******---
%平面机械臂的逆运动学求解，Bruno英文书91页，中文68页,已经本人20180405推导验证
%{
 c2=x*x+y*y-a1*a1-a2*a2; %cos(q2)
 s2=sqrt(1-c2*c2); %对应于q2为正号（肘关节朝上，坐标系中逆时针）
 s2=-sqrt(1-c2*c2); %对应于q2为负号（肘关节朝下）
 q2=atan2(s2,c2);

 s1=((a1+a2*cos(q2))*y-a2*sin(q2)*x)/(x*x+y*y);
 c1=((a1+a2*cos(q2))*x+a2*sin(q2)*y)/(x*x+y*y);
 q1=atan2(s1,c1);
%}
%---****---
% time base vector
  t_f = 4;
  T = (0:Tc:t_f)';
  n = size(T,1); %T矩阵的行数
  %T_move=3; %3秒完成运动1.8米梯形轨迹的运动
  %加速时间0.6s,匀速运动时间[0.6-0.8],
  a1=1; %length
  a2=1; %length
% final time
  x = zeros(n,1);
  y = zeros(n,1);  %水平运动 
  x(1)=0.2;%初始位置
  
  q1 = zeros(n,1);
  q2 = zeros(n,1);
  dq1 = zeros(n,1);
  dq2 = zeros(n,1);
  ddq1 = zeros(n,1);
  ddq2 = zeros(n,1);  
  
  count=2;
  for t=0+Tc:Tc:0.6
	x(count)=x(1) +0.5*(1/0.6)*t*t;
	count=count+1;
  end
  for t=0.6+Tc:Tc:0.8
	x(count)= x(1)+0.5*(1/0.6)*0.6*0.6+1.0*(t-0.6);
	count=count+1;
  end
  
  for t=0.8+Tc:Tc:3
	x(count)= x(1)+0.5*(1/0.6)*0.6*0.6+1.0*(0.8-0.6)+1*(t-0.8)-0.5*(1/(3-0.8))*(t-0.8)^2;
	count=count+1;
  end  
  x(count:n)= x(count-1);
  
  for i=1:n
	[q1(i),q2(i)]=inv_kinematics(x(i),y(i));
  end
  for i=2:n
	dq1(i)=(q1(i)-q1(i-1))/Tc;
	dq2(i)=(q2(i)-q2(i-1))/Tc;
  end
  for i=2:n
	ddq1(i)=(dq1(i)-dq1(i-1))/Tc;
	ddq2(i)=(dq2(i)-dq2(i-1))/Tc;	
  end
  
  % figure(1)
  % plot(dq1,'.'),grid on
  % figure(2)
  % plot(dq2),grid on

  q_d(1:n,:) = [q1 q2];

  dq_d(1:n,:)  = [dq1 dq2];
  ddq_d(1:n,:) = [ddq1 ddq2];

function [q1,q2]=inv_kinematics(x,y)
  a1=1;
  a2=1;
  c2=(x*x+y*y-a1*a1-a2*a2)/(2*a1*a2);  %cos(q2)
  %s2=sqrt(1-c2*c2); %对应于q2为正号（肘关节朝上，坐标系中逆时针）
  s2=sqrt(1-c2*c2); %对应于q2为负号（肘关节朝下）
  q2=atan2(s2,c2);
  s1=((a1+a2*cos(q2))*y-a2*sin(q2)*x)/(x*x+y*y);
  c1=((a1+a2*cos(q2))*x+a2*sin(q2)*y)/(x*x+y*y);
  q1=atan2(s1,c1);
end
%{


%}
%---------end of file---------------