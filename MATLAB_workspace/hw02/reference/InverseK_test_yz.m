model.L1 = 0.5; %length of link1
model.L2 = 0.5; %length of link2//定义杆件长度

w = pi / 2; %旋转角速度
t1 = 2; %t1之前走直线，之后画圈
j = 1;
for t = 0:0.2:8
    pause(0.01)
    if t < t1
        [y, yd, ydd] = TSpline(0, 0, 0, 0, t1/2, 0, 0.25*w, t1, t);
        [z, zd, zdd] = TSpline(-1, 0, 0, -0.9, t1/2, -0.75, 0, t1, t); %三次样条差值，规划末端点的位置、速度、加速度
        % z=-1+0.25*t;
        % zd=0.25;
    else
        y = 0.25 * sin(w*(t - t1));
        yd = 0.25 * w * cos(w*(t - t1));
        ydd = -0.25 * w * w * sin(w*(t - t1));
        z = -0.5 - 0.25 * cos(w*(t - t1));
        zd = 0.25 * w * sin(w*(t - t1));
        zdd = 0.25 * w * w * cos(w*(t - t1)); %用正弦函数规划末端点的位置、速度、加速度
    end
    q2 = acos((y * y + z * z - model.L1 * model.L1 - model.L2 * model.L2)/(2 * model.L1 * model.L2));
    q1 = atan2(y, -z) - acos((model.L1 * model.L1 + y * y + z * z - model.L2 * model.L2)/(2 * model.L1 * sqrt(y*y+z*z))); %逆运动学求关节角度
    % q2=-acos((y*y+z*z-model.L1*model.L1-model.L2*model.L2)/(2*model.L1*model.L2));
    % q1=atan2(y,-z)+acos((model.L1*model.L1+y*y+z*z-model.L2*model.L2)/(2*model.L1*sqrt(y*y+z*z)));%逆运动学的另一组解
    model.Refq(1) = q1;
    model.Refq(2) = q2;
    y1 = model.L1 * sin(q1);
    z1 = -model.L1 * cos(q1); %求L1的末端点的位置
    % yplot=[0 y1 y];
    % zplot=[0 z1 z];
    % axis equal;
    % % axis([-1 1 -1 0])
    % plot(yplot,zplot,'-');
    % hold on
    % plot(y,z,'o')
    % hold on

    plot([0, y1], [0, z1], 'r-')
    hold on
    plot([y1, y], [z1, z], 'b-')
    hold on
    plot(y, z, 'go')
    hold on
    axis equal;
    % plot(y,z,'*')
    % hold off
    % M(j)=getframe;
    % j=j+1;
end
%  movie(M)

% s1=sin(model.Refq(1));
% c1=cos(model.Refq(1));
% s12=sin(model.Refq(1)+model.Refq(2));
% c12=cos(model.Refq(1)+model.Refq(2));
% L1=model.L1;
% L2=model.L2;
% if t<0.003
%     model.Refqd=[0 0]';
%     model.Refqdd=[0 0]';
% else
%     J=[c1*L1+c12*L2 c12*L2;
%        s1*L1+s12*L2 s12*L2];%雅克比矩阵
%     model.Refqd=inv(J)*[yd zd]';%关节角速度
%     model.Refqdd=inv(J)*([ydd
%     zdd]'-[-s1*L1*(model.Refqd(1))^2-s12*L2*(model.Refqd(1)+model.Refqd(2))^2;c1*L1*(model.Refqd(1))^2+c12*L2*(model.Refqd(1)+model.Refqd(2))^2]);;%关节角加速度
% end
