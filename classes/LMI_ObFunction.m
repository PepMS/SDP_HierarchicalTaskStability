classdef (Abstract) LMI_ObFunction
    methods(Abstract)
        [F, bS] = fillLMI(obj);
    end    
end