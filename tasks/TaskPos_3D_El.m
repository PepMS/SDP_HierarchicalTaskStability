classdef TaskPos_3D_El < Task   
    methods
        function J = getTaskJacobian(obj, q_)
            L = obj.robot.links;
            robot = SerialLink(L(1:3), 'name', 'elbow');
            
            J = robot.jacob0(q_(1:3)');
            J([1,3:6], :) = [];
            J = [J, zeros(size(J,1), obj.robot.n - size(J,2))];        
        end
        function e = getTaskError(obj, q_)
            L = obj.robot.links;
            robot = SerialLink(L(1:3), 'name', 'elbow');
            
            r = robot.fkine(q_(1:3));
            e = obj.r_des - r(2,4);
        end
    end
end