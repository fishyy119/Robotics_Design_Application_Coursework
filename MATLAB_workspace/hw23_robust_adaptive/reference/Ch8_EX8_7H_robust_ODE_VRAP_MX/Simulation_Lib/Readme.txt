

//----------------------------------------Note------------------------------------------------//

1.Vrep与MATLAB联合仿真时需要针对自己的电脑进行配置，如不进行该步操作，仿真时会报错！！！
具体配置方法如下：
       三个必备脚本文件：remApi.m、remoteApiProto.m、remoteApi.dll
       依次在Vrep的安装路径下找到这三个必备文件并将之复制、替换Simulation_Lib文件夹里的三个相同文件名的文件。
       其中，remApi.m和remoteApiProto.m两个.m脚本进入Vrep的安装路径：..\CoppeliaRobotics\CoppeliaSimEdu\programming\remoteApiBindings\matlab\matlab路径下寻找；
                 remoteApi.dll文件在..\CoppeliaRobotics\CoppeliaSimEdu\programming\remoteApiBindings\lib\lib\Windows路径下寻找。
       例如，..\CoppeliaRobotics\CoppeliaSimEdu\programming\remoteApiBindings\lib\lib\Windows在我电脑里的完整安装路径就是：C:\Program Files\CoppeliaRobotics\CoppeliaSimEdu\programming\remoteApiBindings\lib\lib\Windows

2. 不要随意改变remApi.m 和 remoteApiProto.m 这两个文件，这是Vrep与Matlab联合仿真的基础文件，系统自带！！！

3.MatlabVrepTemplate.m 包含了整个项目需要调用的所有子函数接口，这是一个可编辑文件。当然，如果你对Vrep仿真尚不熟悉的话，建议不要擅自对其进行更改。

4.twoLinksPlanar.ttt 以及 M1.ttt 文件为项目需要用到的两个仿真模型，分别对应无负载和有负载两种情况。
