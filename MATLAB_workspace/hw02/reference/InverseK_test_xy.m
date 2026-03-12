model.L1 = 0.5; %length of link1
model.L2 = 0.5; %length of link2//定义杆件长度

w = pi / 2; %旋转角速度
t1 = 2; %t1之前走直线，之后画圈
j = 1;
for t = 0:0.2:8
    pause(0.01)
    if t < t1
        [x, xd, xdd] = TSpline(0, 0, 0, 0, t1/2, 0, 0.25*w, t1, t); %[S,Sd,Sdd]=TSpline(p0,v0,t0,p1,t1,p2,v2,t2,t)
        [y, yd, ydd] = TSpline(-1, 0, 0, -0.9, t1/2, -0.75, 0, t1, t); %三次样条差值，规划末端点的位置、速度、加速度
        % z=-1+0.25*t;
        % zd=0.25;
    else
        x = 0.25 * sin(w*(t - t1));
        xd = 0.25 * w * cos(w*(t - t1));
        xdd = -0.25 * w * w * sin(w*(t - t1));
        y = -0.5 - 0.25 * cos(w*(t - t1));
        yd = 0.25 * w * sin(w*(t - t1));
        ydd = 0.25 * w * w * cos(w*(t - t1)); %用正弦函数规划末端点的位置、速度、加速度
    end
    q2 = acos((x * x + y * y - model.L1 * model.L1 - model.L2 * model.L2)/(2 * model.L1 * model.L2));
    q1 = atan2(x, -y) - acos((model.L1 * model.L1 + x * x + y * y - model.L2 * model.L2)/(2 * model.L1 * sqrt(x*x+y*y))); %逆运动学求关节角度
    % q2=-acos((x*x+y*y-model.L1*model.L1-model.L2*model.L2)/(2*model.L1*model.L2));
    % q1=atan2(x,-y)+acos((model.L1*model.L1+x*x+y*y-model.L2*model.L2)/(2*model.L1*sqrt(x*x+y*y)));%逆运动学的另一组解
    model.Refq(1) = q1;
    model.Refq(2) = q2;
    x1 = model.L1 * sin(q1);
    y1 = -model.L1 * cos(q1); %求L1的末端点的位置

    plot([0, x1], [0, y1], 'r-')
    hold on
    plot([x1, x], [y1, y], 'b-')
    hold on
    plot(x, y, 'ro')
    hold on
    axis equal;
    % plot(y,z,'*')
    % hold off
    % M(j)=getframe;
    % j=j+1;
end
%  movie(M)
