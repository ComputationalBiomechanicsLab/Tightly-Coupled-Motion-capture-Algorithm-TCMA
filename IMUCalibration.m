classdef IMUCalibration < AbstractObject
    % IMUCALIBRATION - Generates a calibration (biases, covariances etc.) from
    % given IMU data files.
    
    properties (Description = 'Result')
        GyroscopeBias
        AccelerometerBias
        AccelerometerCalibrationMatrix
        Covariance
    end
    
    properties (Description = 'Files')
        StationaryDataFile
        OrientationDataFile
    end
    
    %% Public methods
    methods
       
        function obj = IMUCalibration(varargin)
            % IMUCalibration constructor
            obj = obj@AbstractObject(varargin{:});
             
        end
        
       function calibrateFromFile(obj)

%% Parse stationary data
allStationaryData = readmatrix(obj.StationaryDataFile);
stationaryData = allStationaryData(:,3:8);
stationaryGyroData = stationaryData(:,4:6);
stationaryAccelData = stationaryData(:,1:3);

% Calculate gyro bias and overall covariance
obj.GyroscopeBias = mean(stationaryGyroData).';
obj.Covariance = cov(stationaryData);

%% Parse orientation data
sensorData = readmatrix(obj.OrientationDataFile);
accelData = sensorData(:,3:5);

[bias,calMatrix] = accelerometerEllipsoidFitting(accelData.');
obj.AccelerometerBias = bias;
obj.AccelerometerCalibrationMatrix = calMatrix;

calibratedAccelData = (calMatrix * (accelData.' - bias)).';

figure();
plot3(accelData(:,1),accelData(:,2),accelData(:,3),'b.'); hold on;
plot3(calibratedAccelData(:,1),calibratedAccelData(:,2),calibratedAccelData(:,3),'k.');
legend('Raw','Calibrated')
xlabel x
ylabel y
zlabel z
axis equal

end
    end
    
end