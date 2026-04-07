function newObj = cut(obj,T)

if ~isempty(obj.RobotData)
    robotData = obj.RobotData;
    robotTime = robotData.Time;
    robotIdx = robotTime > T(1) & robotTime < T(2);
    
    q = robotData.Position;
    qd = robotData.Velocity;
    tau = robotData.Effort;
    
    newRobotData = RobotDataset('Position',q(robotIdx,:),...
        'Velocity',qd(robotIdx,:),'Effort',tau(robotIdx,:),...
        'Time',robotTime(robotIdx));
else
    newRobotData = RobotDataset.empty();
end

if ~isempty(obj.IMUData)
    imuData = obj.IMUData;
    imuTime = imuData(1).Time;
    imuIdx = imuTime > T(1) & imuTime < T(2);
    
    nIMUs = numel(imuData);
    newIMUData(nIMUs) = IMUDataset();
    
    for i = 1:nIMUs
        imu = imuData(i);
        newIMUData(i) = IMUDataset('DeviceID',imu.DeviceID,...
            'Calibration',imu.Calibration,'Time',imu.Time(imuIdx),...
            'CalibratedData',imu.CalibratedData(imuIdx,:));
    end
else
    newIMUData = IMUDataset.empty();
end

if ~isempty(obj.MocapData)
    mocapTime = obj.MocapData.Time;
    mocapIdx = mocapTime > T(1) & mocapTime < T(2);
    
    % Only keep markers
    trackers = obj.MocapData.MotiveData;
    trackers = trackers(strcmpi({trackers.Type},'Marker'));
    
    nMarkers = numel(trackers);
    markerData(nMarkers) = MotiveMarker();
    
    for i = 1:nMarkers
        markerData(i).Time      = trackers(i).Time(mocapIdx);
        markerData(i).Position  = trackers(i).Position(mocapIdx,:);
        markerData(i).Name      = trackers(i).Name;
        markerData(i).Type      = trackers(i).Type;
    end
    
    newMocapData = MotiveDataset('MotiveData',markerData);
else
    newMocapData = MotiveDataset.empty();
end

if ~isempty(obj.ForceData)
    forceData = obj.ForceData;
    forceTime = forceData.Time;
    forceIdx = forceTime > T(1) & forceTime < T(2);
    
    newForceData = ForceDataset('Data',forceData.Data(forceIdx,:),...
        'Time',forceTime(forceIdx),'Covariance',forceData.Covariance);
    
else
    newForceData = ForceDataset.empty();
end

% newName = sprintf('Cut %s',obj.Name);
newName = obj.Name;
newObj = ROSDataset('Name',newName,'MocapData',newMocapData,...
    'IMUData',newIMUData,'RobotData',newRobotData,'ForceData',newForceData);

end