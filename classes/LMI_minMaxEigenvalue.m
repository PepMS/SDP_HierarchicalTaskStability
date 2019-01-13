classdef LMI_minMaxEigenvalue < LMI_ObFunction
    
    methods
        function [F, bS] = fillLMI(obj, A_, F)
            index = min(size(F)) + 1;
            dim = size(A_,1);
            
            bS = dim;
            
            for ii=1:dim
                M_aux = zeros(size(A_,1));
                M_aux(:, ii) = A_(:, ii);
                F{index, ii+1} = -(M_aux + M_aux')/2;
            end
            F{index, ii+2} = eye(size(A_,1));                      
        end
    end
end
