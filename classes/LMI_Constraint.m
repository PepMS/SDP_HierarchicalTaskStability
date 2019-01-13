classdef (Abstract) LMI_Constraint
    methods(Abstract)
        % Constructor
        [F, bS] = fillLMI(obj);
    end    
end