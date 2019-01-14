classdef LMI_jVelBound < LMI_Constraint 
    properties
        u_bound;
        l_bound;    
    end
            
    methods
        %Constructor
        function obj = LMI_jVelBound(u_bound_, l_bound_)            
            obj.u_bound = u_bound_;
            obj.l_bound = l_bound_;
        end
        
        function [F, bS] = fillLMI(obj, Tasks_, F, A_, q_)
            index = size(F,1) + 1;
                             
            dim = length(obj.u_bound);
            bS = [dim, dim];
            
            S = obj.computeS(Tasks_, q_);
            
            for ii=1:size(S,2)
                F{index    , ii+1} =  diag(S(:,ii)); % Lower bound
                F{index + 1, ii+1} = -diag(S(:,ii)); % Upper bound
            end
            F{index, 1} =  diag(obj.l_bound);                      
            F{index + 1, 1} = -diag(obj.u_bound);
        end
        
        function S = computeS(obj, Tasks_, q_)
            S = [];
            
            for ii=1:length(Tasks_)
                task = Tasks_{ii};
                
                % get Jacobian and error
                J = task.getTaskJacobian(q_);
                e = task.getTaskError(q_);
                
                if ii == 1
                    aug_J = J;
                    aug_N = eye(size(J,2));
                else
                    aug_N = eye(size(J,2)) - pinv(aug_J)*aug_J;
                    aug_J = [aug_J; J];
                end
                
                S = [S, aug_N*pinv(J)*diag(e)];
                
            end
        end
    end
end
