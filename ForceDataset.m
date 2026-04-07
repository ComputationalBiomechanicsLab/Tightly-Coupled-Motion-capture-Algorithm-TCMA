classdef ForceDataset < AbstractObject
    % FORCEDATASET - Class that holds data from the 6 DoF force-torque sensor.
    %
    properties (Description = 'Data')
        Data(:,6) double
        Time(:,1) double
    end
    
    properties (Description = 'Calibration')
        Covariance
    end
    
    properties (Dependent, Description = 'Data')
        Force
        Torque
    end
    
    %% Get methods
    methods
        function f = get.Force(obj)
            f = obj.Data(:,1:3);
        end
        
        function M = get.Torque(obj)
            M = obj.Data(:,4:6);
        end
    end
    
    %% Public methods
    methods
        function obj = ForceDataset(varargin)
            % ForceDataset constructor
            obj = obj@AbstractObject(varargin{:});
        end
        
        function [varargout] = plotData(obj,ax)
            if nargin < 2 || isempty(ax)
                ax = linkedSubplots(2,1);
            end
            
            t = obj.Time;
            t = t - t(1);
            f = obj.Force;
            M = obj.Torque;
            
            magf = sqrt(sum(f.^2,2));
            magM = sqrt(sum(M.^2,2));
            
            axes(ax(1));
            h(:,1) = plot(t,f); hold on;
            h(end+1,1) = plot(t,magf,'k--');
            ylabel('Force [N]')
            legend('x','y','z','Magnitude');
            
            axes(ax(2));
            h(1:3,2) = plot(t,M); hold on;
            h(4,2) = plot(t,magM,'k--');
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