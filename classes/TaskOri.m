classdef TaskOri < Task
    methods
        function J = getTaskJacobian(obj, q)
            J = ones(1, length(q));
        end
        function e = getTaskError(obj, q)
            r = sum(q);
            e = obj.r_des - r;
        end
    end
end

