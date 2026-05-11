function motor = getMotorParams(motor_name)
% 返回所选直流电机的 SI 制参数。

switch motor_name
    case 'RE50 200W'
        motor.name = motor_name;
        motor.R = 0.608;
        motor.L = 0.423e-3;
        motor.Kt = 93.4e-3;
        motor.Ke = (1 / 102) * 60 / (2 * pi);
        motor.J = 542e-7;
        motor.U = 48;
    case 'RE40 150W'
        motor.name = motor_name;
        motor.R = 24.4;
        motor.L = 6.41e-3;
        motor.Kt = 266e-3;
        motor.Ke = (1 / 35.9) * 60 / (2 * pi);
        motor.J = 120e-7;
        motor.U = 48;
    case 'RE25 20W'
        motor.name = motor_name;
        motor.R = 12.6;
        motor.L = 1.31e-3;
        motor.Kt = 55e-3;
        motor.Ke = (1 / 174) * 60 / (2 * pi);
        motor.J = 10.5e-7;
        motor.U = 48;
    otherwise
        error('Unsupported motor type: %s', motor_name);
end
end
