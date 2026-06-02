%Example in Section 8.7 
%TextBook: Bruno Siciliano, et al. Robotics: Modelling, Planning and Control, Springer,2009

  a1=1; %arm length
  a2=1; %arm length
  posx_d=a1*cos(q_d(:,1))+a2*cos(q_d(:,1)+q_d(:,2));
  posx=a1*cos(q(:,1))+a2*cos(q(:,1)+q(:,2));
  posy_d=a1*sin(q_d(:,1))+a2*sin(q_d(:,1)+q_d(:,2));
  posy=a1*sin(q(:,1))+a2*sin(q(:,1)+q(:,2));
  norm_error=sqrt((posx_d-posx).^2+(posy_d-posy).^2);
  figure('name','Scheme E: computed torque feedforward');

if c==0,
% joint angle
  fig1=subplot(2,2,1);
  plot(time, q_d(:,1),'--',time,q(:,1)),grid on;
  legend(fig1,{'Ref','Rel'})
  %axis([0 t_d -4000 4000]);
   set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[rad]');
  title('Joint1 angle');

  fig2=subplot(2,2,2);
  plot(time, q_d(:,2),'--',time,q(:,2)),grid on;
  legend(fig2,{'Ref','Rel'})
  %axis([0 t_d -4000 4000]);
   set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[rad]');
  title('Joint2 angle');

% joint torques
 fig3= subplot(2,2,3);
  plot(time, tau(:,1),time,tau(:,2),'--'),grid on;
  legend(fig3,{'Tau1','Tau2'})
  %axis([0 t_d -4000 4000]);
   set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[Nm]');
  title('Joint torque');

% joint position errors
  fig4=subplot(2,2,4);
  %plot(time, (q_d(:,1)-q(:,1)),time, (q_d(:,2)-q(:,2)),'--'),grid on;
  %legend(fig4,{'err1','err2'})
  plot(time, norm_error),grid on;
%  legend(fig4,{'Error_x','Error_y'})
%  axis([0 t_d -0.05 0.05]);
   set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m]');
  title('Norm of position error');


end;