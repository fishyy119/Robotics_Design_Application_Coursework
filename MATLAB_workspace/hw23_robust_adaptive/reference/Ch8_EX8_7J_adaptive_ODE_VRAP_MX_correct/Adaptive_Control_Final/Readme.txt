

//-------------------------Note------------------------//

1. 本文件夹为Adative_Control， 对应自适应控制部分。

2. 仿真模型对应M1.ttt。

3. 其中，params.m, load_m.m, Trapezoid.m ，tra_generation.m 和 Adative_Controller.m为公共文件（或函数）。

4. 执行MAIN_Program.m文件，需要配合Vrep进行仿真，观察控制效果；
    执行ODE_Check.m文件，利用ode113配合DynFunc.m函数进行动力学校验。