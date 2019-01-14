classdef TaskPos_3D_EE < Task
    methods
        function J = getTaskJacobian(obj, q)
            J = obj.robot.jacob0(q');
            J(4:6,:) = [];
        end
        function e = getTaskError(obj, q)
            r = obj.robot.fkine(q);
            % e = obj.r_des - r.t(1:2);
            e = obj.r_des - r(1:3,4);
        end
    end
end

