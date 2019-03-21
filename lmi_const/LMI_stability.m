classdef LMI_stability < LMI_Constraint
    methods
        function [F, bS] = fillLMI(obj, Tasks_, F, A_, q_)
            index = min(size(F)) + 1;
            dim = size(A_,1);
            
            bS = dim;
            
            for ii=1:dim
                A_aux = zeros(size(A_,1));
                A_aux(:, ii) = A_(:, ii);
                F{index, ii+1} = (A_aux + A_aux')/2;
            end
            F{index, 1} =  1*eye(3);
        end
    end
end
