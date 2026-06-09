function traj = build_joint_trajectory(traj_params)
% 使用结构体参数生成关节空间直线梯形速度轨迹。
% traj_params.q_i 为初始关节角列向量。
% traj_params.q_f 为终止关节角列向量。
% traj_params.t_move 为轨迹运动时长。
% traj_params.Ts 为离散采样周期。
% traj_params.t_total 为总仿真时长，可省略，默认等于 t_move。
% traj_params.t_acc_ratio 为加速段占运动时长的比例，可省略，默认 0.25。

q_i = traj_params.q_i;
q_f = traj_params.q_f;
t_move = traj_params.t_move;
Ts = traj_params.Ts;

if isfield(traj_params, 't_total') && ~isempty(traj_params.t_total)
    t_total = traj_params.t_total;
else
    t_total = t_move;
end

if isfield(traj_params, 't_acc_ratio') && ~isempty(traj_params.t_acc_ratio)
    t_acc_ratio = traj_params.t_acc_ratio;
else
    t_acc_ratio = 0.25;
end

if t_acc_ratio <= 0 || t_acc_ratio >= 0.5
    error('t_acc_ratio must satisfy 0 < t_acc_ratio < 0.5.');
end

time_move = (0:Ts:t_move)';
t_acc = t_acc_ratio * t_move;
n_move = numel(time_move);

s = zeros(n_move, 1);
ds = zeros(n_move, 1);
dds = zeros(n_move, 1);

v_c = 1 / (t_move - t_acc);
a_c = v_c / t_acc;

for i = 1:n_move
    ti = time_move(i);
    if ti <= t_acc
        s(i) = 0.5 * a_c * ti^2;
        ds(i) = a_c * ti;
        dds(i) = a_c;
    elseif ti <= t_move - t_acc
        s(i) = 0.5 * a_c * t_acc^2 + v_c * (ti - t_acc);
        ds(i) = v_c;
    else
        dt_dec = t_move - ti;
        s(i) = 1 - 0.5 * a_c * dt_dec^2;
        ds(i) = a_c * dt_dec;
        dds(i) = -a_c;
    end
end

s(1) = 0;
s(end) = 1;
ds(1) = 0;
ds(end) = 0;
dds([1, end]) = 0;

delta_q = q_f - q_i;

q_move = q_i + delta_q * s';
dq_move = delta_q * ds';
ddq_move = delta_q * dds';

% 补全末尾常值段
if t_total <= t_move
    traj.t = time_move;
    traj.q = q_move;
    traj.dq = dq_move;
    traj.ddq = ddq_move;
    return;
end

time_total = (0:Ts:t_total)';
n_total = numel(time_total);

traj.t = time_total;
traj.q = [q_move, repmat(q_f, 1, n_total-n_move)];
traj.dq = [dq_move, zeros(numel(q_i), n_total-n_move)];
traj.ddq = [ddq_move, zeros(numel(q_i), n_total-n_move)];
end
