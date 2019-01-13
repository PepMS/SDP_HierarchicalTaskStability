classdef TaskPos < Task
    methods
        function J = getTaskJacobian(obj, q)
            J = obj.robot.jacob0(q');
            J(3:6,:) = [];
        end
        function e = getTaskError(obj, q)
            r = obj.robot.fkine(q);
            e = obj.r_des - r.t(1:2);
            % e = obj.r_des - r(1:2,4);
        end
    end
end

