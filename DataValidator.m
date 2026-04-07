classdef DataValidator < AbstractObject
    
    properties (Description = 'Input')
        MeasurementData(1,1) ROSDataset
        EstimatedData(1,1) RobotDataset
    end
    
    properties (Dependent,Description = 'Input')
        Q
        QDot
        Tau
    end
    
    properties (Transient, Description = 'Model')
        Model(1,1) OpenSimModel
    end
    
    properties (Description = 'Model')
        SensorMapping
    end
    
    properties (Description = 'Output')
        ModelMarkerData
        ModelIMUData
        MarkerError
        IMUError
    end
    
    properties (Dependent, Description = 'Output')
        AngularVelocityError
        LinearAccelerationError
        QError
        QDotError
        TauError
        QRMSE
        QDotRMSE
        TauRMSE
        MarkerRMSE
        AngularVelocityRMSE
        LinearAccelerationRMSE
    end
    
    %% Get methods
    methods
        function q = get.Q(obj)
            q = obj.EstimatedData.Position;
        end
        
        function qd = get.QDot(obj)
            qd = obj.EstimatedData.Velocity;
        end
        
        function tau = get.Tau(obj)
            tau = obj.EstimatedData.Effort;
        end
        
        function we = get.AngularVelocityError(obj)
            if ~isempty(obj.IMUError)
                we = [obj.IMUError.AngularVelocity];
            else
                we = [];
            end
        end
        
        function ae = get.LinearAccelerationError(obj)
            if ~isempty(obj.IMUError) 
                ae = [obj.IMUError.LinearAcceleration];
            else
                ae = [];
            end
        end
        
        function qe = get.QError(obj)
            idxFree = obj.Model.idxFree;
            if size(obj.Q,2) > numel(obj.Model.FreeCoordinates)
                q = obj.Q(:,idxFree);
            else
                q = obj.Q;
            end
            if ~isempty(idxFree) && ~isempty(obj.MeasurementData.RobotData) && ~isempty(obj.Q)
                qe = obj.MeasurementData.RobotData.Position(:,idxFree) - q;
            else
                qe = [];
            end
        end
        
        function qde = get.QDotError(obj)
            idxFree = obj.Model.idxFree;
            if size(obj.QDot,2) > numel(obj.Model.FreeCoordinates)
                qd = obj.QDot(:,idxFree);
            else
                qd = obj.QDot;
            end
            if ~isempty(idxFree) && ~isempty(obj.MeasurementData.RobotData) && ~isempty(obj.QDot)
                qde = obj.MeasurementData.RobotData.Velocity(:,idxFree) - qd;
            else
                qde = [];
            end
        end
        
        function taue = get.TauError(obj)
            idxFree = obj.Model.idxFree;
            if size(obj.Tau,2) > numel(obj.Model.FreeCoordinates)
                tau = obj.Tau(:,idxFree);
            else
                tau = obj.Tau;
            end
            if ~isempty(idxFree) && ~isempty(obj.MeasurementData.RobotData) && ~isempty(obj.Tau)
                taue = -obj.MeasurementData.RobotData.Effort(:,idxFree) - tau;
            else
                taue = [];
            end
        end
        
        function qe = get.QRMSE(obj)
            qe = rms(obj.QError);
        end
        
        function qde = get.QDotRMSE(obj)
            qde = rms(obj.QDotError);
        end
        
        function taue = get.TauRMSE(obj)
            taue = rms(obj.TauError);
        end
        
        function me = get.MarkerRMSE(obj)
            if ~isempty(obj.MarkerError)
                nMarkers = numel(obj.MarkerError);
                me = zeros(nMarkers,1);
                for i = 1:nMarkers
                    % Take rms of norm
                    d = sqrt(sum(obj.MarkerError(i).Position.^2,2));
                    me(i) = rms(d,'omitnan');
                end
            else
                me = [];
            end
        end
        
        function we = get.AngularVelocityRMSE(obj)
            if ~isempty(obj.AngularVelocityError)
                we = zeros(numel(obj.IMUError),1);
                for i = 1:numel(obj.IMUError)
                    idx = (i-1)*3 + (1:3);
                    e = obj.AngularVelocityError(:,idx);
                    d = sqrt(sum(e.^2,2));
                    we(i) = rms(d,'omitnan');
                end
            else
                we = [];
            end
        end
        
        function ae = get.LinearAccelerationRMSE(obj)
            if ~isempty(obj.LinearAccelerationError)
                ae = zeros(numel(obj.IMUError),1);
                for i = 1:numel(obj.IMUError)
                    idx = (i-1)*3 + (1:3);
                    e = obj.LinearAccelerationError(:,idx);
                    d = sqrt(sum(e.^2,2));
                    ae(i) = rms(d,'omitnan');
                end
            else
                ae = [];
            end
        end
    end
    
    %% Constructor
    methods
        function obj = DataValidator(varargin)
            obj = obj@AbstractObject(varargin{:});
        end
    end
    
    %% Public methods
    methods
       function calculateMeasurement(obj)
import org.opensim.modeling.*

model = obj.Model;
model.osObj.initSystem();

markerData = obj.MeasurementData.MocapData.MotiveData;

if ~isempty(markerData)
    markerData = markerData.select({model.Markers.Name});
end

imuData = obj.MeasurementData.IMUData;
robotData = obj.MeasurementData.RobotData;

t = robotData.Time;
nt = numel(t);
dt = mean(diff(t));

nMarkers = numel(markerData);
nIMUs = numel(imuData);

Y_marker = zeros(nt,3*nMarkers);
Y_imu = zeros(nt,6*nIMUs);

if size(obj.Q,2) > numel(model.FreeCoordinates)
    idx = model.idxFree;
    Q = obj.Q(:,idx);
    Qd = obj.QDot(:,idx);
    Tau = obj.Tau(:,idx);
else
    Q = obj.Q;
    Qd = obj.QDot;
    Tau = obj.Tau;
end

% Qdd = finiteDifference(Qd,dt);
Qdd = zeros(size(Qd));

for i = 1:nt
    qi = Q(i,:);
    qdi = Qd(i,:);
    % qddi = Qdd(i,:);
    taui = Tau(i,:);
    ti = t(i);
    
    qddi = model.calcQdd(qi,qdi,taui,ti);
    Qdd(i,:) = qddi ;
    
    Y_marker(i,:) = model.getMarkerPositions(qi);
    Y_imu(i,:) = model.getIMUMeasurement(qi(:),qdi(:),qddi(:));
end

modelMarkerData(nMarkers) = MotiveMarker();
markerError(nMarkers) = MotiveMarker();

for i = 1:nMarkers
    idx = (i-1)*3 + (1:3);
    modelMarkerData(i).Name = markerData(i).Name;
    modelMarkerData(i).Time = t;
    modelMarkerData(i).Type = markerData(i).Type;
    
    modelMarkerData(i).Position = Y_marker(:,idx);
    
    markerError(i).Name = markerData(i).Name;
    markerError(i).Time = t;
    markerError(i).Type = markerData(i).Type;
    
    markerError(i).Position = Y_marker(:,idx) - markerData(i).Position;
end

modelIMUData(nIMUs) = IMUDataset();
imuError(nIMUs) = IMUDataset();

for i = 1:nIMUs
    idx = (i-1)*6 + (1:6);
    modelIMUData(i).DeviceID = model.IMUs(i).Name;
    modelIMUData(i).Time = t;
    
    modelIMUData(i).CalibratedData = Y_imu(:,idx);
    
    imuError(i).DeviceID = model.IMUs(i).Name;
    imuError(i).Time = t;
    
    deviceIdx = strcmp(obj.SensorMapping.('IMU Name'),model.IMUs(i).Name);
    deviceID = obj.SensorMapping{deviceIdx,'Device ID'};
    y_meas = imuData.getData(deviceID);
    imuError(i).CalibratedData = Y_imu(:,idx) - y_meas;
end

obj.ModelMarkerData = modelMarkerData;
obj.ModelIMUData = modelIMUData;

obj.MarkerError = markerError;
obj.IMUError = imuError;

end
    end
    
end