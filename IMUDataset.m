classdef IMUDataset < AbstractObject
    % IMUDATASET - Class for storing IMU data (i.e. linear acceleration and
    % angular velocity), as well as calibration parameters etc.
    
    properties (Description = 'Data')
        RawData(:,6) double
        CalibratedData(:,6) double
        Time(:,1) double
    end
    
%     properties (Description = 'Calibration')
%         AccelerationCalibration(3,3) double = eye(3)
%         AccelerationBias(3,1) double
%         GyroscopeBias(3,1) double
%         Covariance(6,6) double
%     end
    
    properties (Description = 'Device')
        DeviceID {char,string}
        Calibration(1,1) IMUCalibration
    end
    
    properties (Dependent, Description = 'Data')
        LinearAcceleration
        AngularVelocity
    end
    
    %% Get methods
    methods
        function a = get.LinearAcceleration(obj)
            a = obj.CalibratedData(:,1:3);
        end
        
        function w = get.AngularVelocity(obj)
            w = obj.CalibratedData(:,4:6);
        end
    end
    
    %% Public methods
    methods
        function objsel = select(objs,name)
            names = {objs.DeviceID};
            selected = contains(names,name);
            objsel = objs(selected);
        end
        
        function obj = IMUDataset(varargin)
            % IMUDataset constructor
            obj = obj@AbstractObject(varargin{:});
        end
        
        function calibrateData(obj)
            
            % Recursive call for arrays
            if numel(obj) > 1
                for i = 1:numel(obj)
                    obj(i).calibrateData();
                end
                return
            end
            
            aRaw = obj.RawData(:,1:3);
            wRaw = obj.RawData(:,4:6);
            
            bAccel = obj.Calibration.AccelerometerBias;
            D = obj.Calibration.AccelerometerCalibrationMatrix;
            bGyro = obj.Calibration.GyroscopeBias;
            
            aCal = (D * (aRaw.' - bAccel)).';
            wCal = wRaw - bGyro.';
            
            obj.CalibratedData = [aCal, wCal];
        end
        
        function [varargout] = plotData(obj,ax)
            if nargin < 2 || isempty(ax)
                ax = linkedSubplots(2,1);
            end
            
            t = obj.Time;
            t = t - t(1);
            a = obj.LinearAcceleration;
            maga = sqrt(sum(a.^2,2));
            
            w = obj.AngularVelocity;
            magw = sqrt(sum(w.^2,2));
            
            axes(ax(1));
            h(:,1) = plot(t,a); hold on;
            h(end+1,1) = plot(t,maga,'k--');
            ylabel('Specific force [m/s$^2$]')
            legend('x','y','z','Magnitude');
            
            axes(ax(2));
            h(1:3,2) = plot(t,w); hold on;
            h(4,2) = plot(t,magw,'k--');
            ylabel('Angular velocity [rad/s]')
            xlabel('Time [s]');
            
            h = [h(:)];
            
            if nargout >= 1
                varargout{1} = h;
            end
            if nargout >= 2
                varargout{2} = ax;
            end
        end
        
        function Y = getData(objs,idx)
            if nargin > 1 && ~isempty(idx)
                if iscell(idx) || ischar(idx) || isstring(idx)
                    names = {objs.DeviceID};
                    [~,idx] = ismember(idx,names);
                end
                objs = objs(idx);
            end
            
            Y = [objs.CalibratedData];
        end
    end
    
end