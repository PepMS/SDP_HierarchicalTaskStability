classdef LMI_gainBound < LMI_Constraint 
    methods
        function [F, bS] = fillLMI(obj, Tasks_, F, M_)
            index = size(F,1) + 1;
            
            [l_bound, u_bound] = obj.getBounds(Tasks_);
                        
            dim = length(u_bound);
            bS = [dim, dim];
            
            for ii=1:dim
                M_aux = zeros(dim);
                M_aux(ii, ii) = 1;
                F{index    , ii+1} =  M_aux; % Lower bound
                F{index + 1, ii+1} = -M_aux; % Upper bound
            end
            F{index, 1} =  diag(l_bound);                      
            F{index + 1, 1} = -diag(u_bound);
        end
        
        function [l_b, u_b] = getBounds(obj, Tasks_)
            l_b = [];
            u_b = [];
            for ii=1:size(Tasks_,2)
                task = Tasks_{ii};
                l_b = [l_b; task.bound_l];
                u_b = [u_b; task.bound_u];
            end            
        end
        
        
    end
end
