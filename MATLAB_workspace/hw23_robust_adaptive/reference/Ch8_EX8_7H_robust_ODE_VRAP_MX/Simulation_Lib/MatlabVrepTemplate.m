classdef MatlabVrepTemplate
    
    properties
        Joint_Num;
        Joint_Name;
        Joint_Handle;
        Body_Handle;
        Real_JointAgl;
        Joint_Agl;
        rJoint;
        Port;
        ClientID;
        Control_Time;
        Main;
    end
    
    methods
        function mod = MatlabVrepTemplate(T_control)
            mod.Joint_Num = 2;
            mod.Control_Time = T_control;
            mod.Joint_Name = cell(1, mod.Joint_Num);
            mod.Joint_Name = {              
                {'Joint_1'},{'Joint_2'}
                };

            mod.Joint_Handle = zeros(1, mod.Joint_Num);
            mod.Real_JointAgl = zeros(1,mod.Joint_Num);
            mod.Joint_Agl = zeros(1, mod.Joint_Num);
            mod.Port = 19997;
            mod.Main = remApi('remoteApi');
        end
          
        function mod2 = init(mod)
            mod2 = mod;
            mod2.Main.simxFinish(-1);
            mod2.ClientID = mod2.Main.simxStart('127.0.0.1',mod2.Port,true,true,5000,mod2.Control_Time);
            for i = 1:mod2.Joint_Num
                Jointname = mod2.Joint_Name{i};
                if i > 0
                    Jointname = Jointname{1};
                end
                [~, mod2.Joint_Handle(i)] = mod2.Main.simxGetObjectHandle(mod2.ClientID,Jointname,mod2.Main.simx_opmode_oneshot_wait);
            end
            [~, mod2.Body_Handle] = mod2.Main.simxGetObjectHandle(mod2.ClientID,'Base',mod2.Main.simx_opmode_blocking);
           end
        
        function go(mod)
            mod.Main.simxSynchronous(mod.ClientID,true);
            mod.Main.simxStartSimulation(mod.ClientID,mod.Main.simx_opmode_blocking);
            mod.trigger();
            mod.trigger();
        end
        
        function trigger(mod)
            mod.Main.simxSynchronousTrigger(mod.ClientID);
        end
        
        function stop(mod)
            mod.Main.simxStopSimulation(mod.ClientID,mod.Main.simx_opmode_blocking);
            disp('Simulation stopped');
        end
        
        function set_joint(mod)
            for i = 1:mod.Joint_Num
                mod.Main.simxSetJointTargetPosition(mod.ClientID,mod.Joint_Handle(i),mod.Joint_Agl(i),mod.Main.simx_opmode_oneshot);
            end
        end
        
        function [rJoint] = read_joint(mod)
            rJoint = mod.Real_JointAgl;
            for i=1:mod.Joint_Num
                [~, rJoint(i)] = mod.Main.simxGetJointPosition(mod.ClientID,mod.Joint_Handle(i),mod.Main.simx_opmode_oneshot);
            end
        end
        
        function enable_joint_torque(mod)
            for i = 1:mod.Joint_Num
                mod.Main.simxSetObjectIntParameter(mod.ClientID, mod.Joint_Handle(i),int32(2001),0,mod.Main.simx_opmode_oneshot);                               
            end            
        end
        
        function set_joint_maxTorque(mod)
            for i = 1:mod.Joint_Num
                mod.Main.simxSetJointForce(mod.ClientID, mod.Joint_Handle(i), 100000, mod.Main.simx_opmode_oneshot);                
            end            
        end
            
    end
end
            