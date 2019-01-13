classdef SDPCLIKProblem
    properties
        robot;   % Serial Link manipulator
        q0;      % Initial robot configuration
        Tasks;   % Array with all tasks to be solved
        dt;      % Delta time I want to consider (in s)
        t_end;   % Ending time of the simulation
        dim_err; % Dimension of the augmented error vector
        
        OF_LMI;  % Object that detrmines the LMI linked with the OF
        LMI_l;   % List of LMIs that has to be added as constraints
        
        AA; % Matrix A container
        MM; % Matrix M container
        eAA; % Matrix A eigenvalues container
        eMM; % Matrix M eigenvalues container
        
    end
    
    methods
        function obj = SDPCLIKProblem(robot_, q0_, Tasks_, dt_, t_end_)
            obj.robot = robot_;
            obj.q0 = q0_;
            obj.Tasks = Tasks_;
            obj.dt = dt_;
            obj.t_end = t_end_;
            obj.dim_err = obj.computeAugmentedDimension();
        end
        
        function dim = computeAugmentedDimension(obj)
            dim = 0;
            for ii = 1:length(obj.Tasks)
                task = obj.Tasks{ii};
                dim = dim + task.dim;
            end
        end
        
        function solve(obj)
            tt = 0:obj.dt:obj.t_end;
            q = obj.q0;
            for t=tt
                % Compute Gains
                K = obj.computeGains(q);
                
                % Compute joint velocity & position
                q_d = obj.computeJVel(K, q);
                q = q + q_d*obj.dt;
                
            end
        end
        
        function K = computeGains(obj, q_)
            
            A = obj.computeA(q_);
            
            nVars = 1 + obj.dim_err;
            blockStruct = [];
            F = {};
            
            % O.F. LMI
            [F, dim] = obj.OF_LMI.fillLMI(A, F);
            blockStruct = [blockStruct, dim];
            
            % Pure LMIs
            for ii=1:length(obj.LMI_l)
                lmi = obj.LMI_l{ii};
                [F, dim] = lmi.fillLMI(obj.Tasks, F, A);
                blockStruct = [blockStruct, dim];
            end
            
            
            c = zeros(1, nVars);
            c(end) = 1;
            
            nBlock = size(F,1);
            OPTION.print = '';
            [objVal, xOpt, X, Y, INFO] = sdpam(nVars, nBlock, blockStruct, c, F, OPTION);
            K = double(xOpt(1:end-1));
        end
        
        function q_d = computeJVel(obj, K_, q_)
            dim_a = 1;
            q_d = zeros(size(q_));
            
            for ii=1:length(obj.Tasks)
                task = obj.Tasks{ii};
                
                % Get gains corresponding to this task
                k = K_(dim_a:dim_a + task.dim - 1);
                dim_a = dim_a + task.dim;
                
                % get Jacobian and error
                J = task.getTaskJacobian(q_);
                e = task.getTaskError(q_);
                
                if ii == 1
                    aug_J = J;
                    aug_N = eye(obj.robot.n);
                else
                    aug_N = eye(obj.robot.n) - pinv(aug_J)*aug_J;
                    aug_J = [aug_J; J];
                end
                
                q_d = q_d + aug_N*pinv(J)*diag(k)*e;
            end
            
        end
        
        function A = computeA(obj, q_)
            A = zeros(obj.dim_err);
            col = 1;
            for ii = 1:length(obj.Tasks)
                post_Task = obj.Tasks{ii};
                post_J  = post_Task.getTaskJacobian(q_);
                
                if ii == 1
                    aug_J = post_J;
                    aug_N = eye(obj.robot.n);
                else
                    aug_N = eye(obj.robot.n) - pinv(aug_J)*aug_J;
                    aug_J = [aug_J; post_J];
                end
                
                row = 1;
                for jj = 1:length(obj.Tasks)
                    pre_Task = obj.Tasks{jj};
                    pre_J = pre_Task.getTaskJacobian(q_);
                    
                    A(row:row+pre_Task.dim-1, col:col+post_Task.dim-1) = pre_J*aug_N*pinv(post_J);
                    
                    row = row + pre_Task.dim;
                end
                col = col + post_Task.dim;
            end
        end
        
    end
end

