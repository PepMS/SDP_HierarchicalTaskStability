classdef LMI_stability < LMI_Constraint
    methods
        function [F, bS] = fillLMI(obj, Tasks_, F, M_)
            index = min(size(F)) + 1;
            dim = size(M_,1);
            
            bS = dim;
            
            for ii=1:dim
                M_aux = zeros(size(M_,1));
                M_aux(:, ii) = M_(:, ii);
                F{index, ii+1} = (M_aux + M_aux')/2;
            end
            F{index, 1} =  eye(3);
        end
    end
end
