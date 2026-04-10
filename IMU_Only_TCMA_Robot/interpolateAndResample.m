function newObj = interpolateAndResample(obj,fs)

if nargin < 2 || isempty(fs)
    fs = 100;
end

if ~isempty(obj.RobotData)
    robotTime = obj.RobotData.Time;
else
    robotTime = [];
end

if ~isempty(obj.IMUData)
    imuTime = obj.IMUData(1).Time;
else
    imuTime = [];
end

if ~isempty(obj.MocapData)
    mocapTime = obj.MocapData.Time;
else
    mocapTime = [];
end

if ~isempty(obj.ForceData)
    forceTime = obj.ForceData.Time;
else
    forceTime = [];
end

% Only take data where all measurements are available
t0 = max([min(robotTime),min(imuTime),min(mocapTime),min(forceTime)]);
t_end = min([max(robotTime),max(imuTime),max(mocapTime),max(forceTime)]);

dt = 1/fs;
t_new = t0:dt:t_end;

nIMUs = numel(obj.IMUData);

imuData(nIMUs) = IMUDataset();

for i = 1:nIMUs
    data = obj.IMUData(i).CalibratedData;
    newData = interp1(imuTime,data,t_new);
    
    imuData(i).Time = t_new;
    imuData(i).CalibratedData = newData;
    imuData(i).Calibration = obj.IMUData(i).Calibration;
    imuData(i).DeviceID = obj.IMUData(i).DeviceID;
end

% Only keep markers, since interpolating rotations is hard
trackers = obj.MocapData.MotiveData;
trackers = trackers(strcmpi({trackers.Type},'Marker'));

nMarkers = numel(trackers);
markerData(nMarkers) = MotiveMarker();

for i = 1:nMarkers
    data = trackers(i).Position;
    newData = interp1(mocapTime,data,t_new);
    
    markerData(i).Time = t_new;
    markerData(i).Position = newData;
    markerData(i).Name = trackers(i).Name;
    markerData(i).Type = trackers(i).Type;
end

newMocapData = MotiveDataset('MotiveData',markerData);
newName = sprintf('%s@%dHz',obj.Name,fs);
newObj = ROSDataset('Name',newName,'MocapData',newMocapData,'IMUData',imuData);

if ~isempty(obj.RobotData)
    q = obj.RobotData.Position;
    qd = obj.RobotData.Velocity;
    tau = obj.RobotData.Effort;
    newQ = interp1(robotTime,q,t_new);
    newQd = interp1(robotTime,qd,t_new);
    newTau = interp1(robotTime,tau,t_new);
    
    robotData = RobotDataset('Position',newQ,'Velocity',newQd,'Effort',newTau,'Time',t_new);
    newObj.RobotData = robotData;
end

if ~isempty(obj.ForceData)
    data = obj.ForceData.Data;
    newData = interp1(forceTime,data,t_new);
    
    forceData = ForceDataset('Data',newData,'Time',t_new,'Covariance',...
        obj.ForceData.Covariance);
    newObj.ForceData = forceData;
end

end