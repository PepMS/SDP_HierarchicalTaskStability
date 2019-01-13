classdef LMI_minMaxEigenvalue < LMI_ObFunction
    
    methods
        function F = fillLMI(obj, M_, F)
            index = min(size(F)) + 1;
            dim = size(M_,1);
            
            for ii=1:dim
                M_aux = zeros(size(M_,1));
                M_aux(:, ii) = M_(:, ii);
                F{index, ii+1} = M_aux + M_aux' - diag(diag(M_aux));
            end
            F{index, ii+2} = eye(size(M_,1));                      
        end
    end
end
