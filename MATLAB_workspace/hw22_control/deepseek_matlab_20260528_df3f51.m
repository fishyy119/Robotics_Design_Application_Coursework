%% 自动生成一个使用 MATLAB Function 模块的 Simulink 模型 (robot_dynamics_example.mdl)
% 模型功能：输入 tau_diff (2x1) 和 cos_q2 (标量)，输出关节加速度 ddq (2x1)
% 需要 BMatrix.m 在 MATLAB 路径上

% 创建新模型
modelName = 'robot_dynamics_example';
new_system(modelName);
open_system(modelName);

% 添加模块
% 1. 两个 Constant 作为输入源
tau_diff_val = [1; 0.5];   % 示例值
cos_q2_val = 0.5;

add_block('simulink/Sources/Constant', [modelName '/tau_diff'], 'Value', mat2str(tau_diff_val));
add_block('simulink/Sources/Constant', [modelName '/cos_q2'], 'Value', num2str(cos_q2_val));

% 2. MATLAB Function 模块
add_block('simulink/User-Defined Functions/MATLAB Function', [modelName '/Compute Acceleration']);

% 3. Scope 显示输出
add_block('simulink/Sinks/Scope', [modelName '/ddq']);

% 设置 MATLAB Function 模块的代码
% 获取模块句柄并设置函数脚本
hFunc = get_param([modelName '/Compute Acceleration'], 'handle');
% 注意：设置 MATLAB Function 的代码需要使用 set_param('...', 'Script', '...')
% 这里提供函数内容
funcScript = [
    'function ddq = fcn(tau_diff, cos_q2)\n', ...
    '    % 计算关节加速度: ddq = inv(B(q)) * tau_diff\n', ...
    '    B = BMatrix(cos_q2);  % 调用已有函数计算惯性矩阵\n', ...
    '    ddq = B \\ tau_diff;   % 左除，效率更高且数值稳定\n', ...
    'end'
];
set_param([modelName '/Compute Acceleration'], 'Script', funcScript);

% 设置输入输出端口的数据属性 (在 MATLAB Function 模块中)
% 通过 set_param 设置端口维度和类型，较复杂，推荐手动修改，但为了完整，使用以下命令
% 注意：这些命令可能因版本略有不同，运行后若出错请手动调整端口维度。
try
    % 设置输入端口 1 (tau_diff) 为 2x1 双精度
    set_param([modelName '/Compute Acceleration'], 'Inputs', 'tau_diff, cos_q2');
    set_param([modelName '/Compute Acceleration'], 'Outputs', 'ddq');
    % 设置端口维度：需要用到 PortDimensions 参数，但不同版本语法有差异，这里跳过自动设置
    % 用户可能需要在打开模型后手动检查端口设置（双击模块 -> Edit Data）。
    disp('模型创建成功。请双击 MATLAB Function 模块，检查端口维度和数据类型：');
    disp(' - 输入 tau_diff: size [2,1], type double');
    disp(' - 输入 cos_q2:   size 1,    type double');
    disp(' - 输出 ddq:      size [2,1], type double');
catch
    disp('请注意：由于 MATLAB 版本差异，可能需要手动配置模块端口属性（双击模块 -> Edit Data）。');
end

% 连接模块
add_line(modelName, 'tau_diff/1', 'Compute Acceleration/1');
add_line(modelName, 'cos_q2/1', 'Compute Acceleration/2');
add_line(modelName, 'Compute Acceleration/1', 'ddq/1');

% 保存模型为 .mdl 格式
save_system(modelName, [modelName '.mdl']);

% 打开模型供查看
open_system(modelName);

disp(['模型已保存为 ' modelName '.mdl，请确保 BMatrix.m 在 MATLAB 路径中。']);