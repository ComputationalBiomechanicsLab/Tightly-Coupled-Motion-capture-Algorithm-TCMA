classdef RobotDataset < AbstractObject
    % ROBOTDATASET - Class that holds joint encoder & torque measurements
    % from a robot
    %
    properties (Description = 'Data')
        Position
        Velocity
        Effort
        Time
    end
    
    methods
        function obj = RobotDataset(varargin)
            % RobotDataset constructor
            obj = obj@AbstractObject(varargin{:});
        end
        
        function P = calcPower(obj)
            P = obj.Effort .* obj.Velocity;
        end
        
        function [varargout] = plotData(obj,ax,style,idx,T)
            if nargin < 2 || isempty(ax)
                ax = linkedSubplots(3,1);
            end
            if nargin < 3 || isempty(style)
                style = '-';
            end
            if nargin < 4 || isempty(idx)
                idx = ':';
            end
            if nargin < 5 || isempty(T)
                T = [-inf inf];
            end
            
            t = obj.Time;
            t = t - t(1);
            
            idx_t = t > T(1) & t < T(2);
            t = t(idx_t);
            
            q = rad2deg(obj.Position(idx_t,idx));
            qd = rad2deg(obj.Velocity(idx_t,idx));
            tau = obj.Effort(idx_t,idx);
            
            axes(ax(1));
            h(:,1) = plot(t,q,style);
            ylabel('Angle [deg]')
            
            axes(ax(2));
            h(:,2) = plot(t,qd,style);
            ylabel('Velocity [deg/s]')
            
            axes(ax(3));
            h(:,3) = plot(t,tau,style);
            ylabel('Torque [Nm]')
            xlabel('Time [s]');
            
            h = [h(:)];
            
            if nargout >= 1
                varargout{1} = h;
            end
            if nargout >= 2
                varargout{2} = ax;
            end
        end
    end
end