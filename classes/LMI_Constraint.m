classdef (Abstract) LMI_Constraint
    methods(Abstract)
        % Constructor
        F = fillLMI(obj);
    end    
end