classdef TaskPos_3D_Wr < Task   
    methods
        function J = getTaskJacobian(obj, q_)
            L = obj.robot.links;
            robot = SerialLink(L(1:4), 'name', 'wrist');
            
            J = robot.jacob0(q_(1:4)');
            J([1,3:6], :) = [];
            J = [J, zeros(size(J,1), obj.robot.n - size(J,2))];
        end
        function e = getTaskError(obj, q_)
            L = obj.robot.links;
            robot = SerialLink(L(1:4), 'name', 'wrist');
            
            r = robot.fkine(q_(1:4));
            % e = obj.r_des - r.t(1:2);
            e = obj.r_des - r(2,4);
        end
    end
end