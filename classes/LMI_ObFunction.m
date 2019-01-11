classdef (Abstract) LMI_ObFunction
    methods(Abstract)
        F = fillLMI(obj);
    end    
end