classdef (Abstract) Task
    properties
        robot; % Task Owner
        dim;   % Task space dimension
        pri;   % Priority of the task
        r_des; % Desired task value
        bound_u; % Gain Upper bound
        bound_l; % Gain Lower bound
    end
    
    methods
        % Constructor
        function obj = Task(robot_, pri_, r_des_)
            try
                if isa(robot_, 'SerialLink')
                    obj.robot = robot_;
                else
                    error(' ');
                end
            catch
                error('Pass, at least, a SerialLink object as a first argument');
            end
            if nargin == 1
                obj.pri = 1;
                obj.dim = 1;
                obj.r_des = zeros(obj.dim, 1);
            elseif nargin == 2
                obj.pri = pri_;
                obj.dim = 1;
                obj.r_des = zeros(obj.dim, 1);
            else
                obj.dim = length(r_des_);
                obj.pri = pri_;
                obj.r_des = r_des_;
            end
        end
    end
    
    % Virtual methods
    methods(Abstract)
        J = getTaskJacobian(obj, q);
        e = getTaskError(obj, q);
    end
end
