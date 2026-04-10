function dataset = readROSBag(file)

%% Get ROS bag object
rosBagSelect = rosbag(file);
availableTopics = rosBagSelect.AvailableTopics.Row;
sec0 = [];

[~,fn,~] = fileparts(file);
dataset = ROSDataset('Name',fn);

%% Read IMUs
% imuData = struct();

if ismember('/xsens_node_01/data',availableTopics)
    imuSelect = rosBagSelect.select("Topic",'xsens_node_01/data');
    imuMsgs = readMessages(imuSelect,'DataFormat','struct');
    imuMsgs = [imuMsgs{:}];
    
    nMsgs = numel(imuMsgs);
    
    % Break apart first message
    msgArray = imuMsgs(1).Mimus;
    nIMUs = numel(msgArray);
    
    imuTime = zeros(nMsgs,1);
    dtIMU = imuMsgs(2).Mimus(1).Imu.Header.Stamp.Nsec - imuMsgs(1).Mimus(1).Imu.Header.Stamp.Nsec;
    dtIMU = double(dtIMU)/1e9;
    
    if isempty(sec0)
        sec0 = msgArray(1).Imu.Header.Stamp.Sec;
    end
    
    imuData(nIMUs) = IMUDataset;
    
    for j = 1:nIMUs
        msg = msgArray(j).Imu;
        deviceID = msg.Header.FrameId;
        imuData(j).DeviceID = deviceID;
        
        % Initialize arrays
        imuData(j).RawData = zeros(nMsgs,6);
        
        acc = [msg.LinearAcceleration.X, msg.LinearAcceleration.Y, msg.LinearAcceleration.Z];
        omega = [msg.AngularVelocity.X, msg.AngularVelocity.Y, msg.AngularVelocity.Z];
        
        % Fill first element
        imuData(j).RawData(1,:) = [acc, omega];
        
        % Set time array & fill first element
        imuData(j).Time = imuTime;
        nSecs = msgArray(j).Imu.Header.Stamp.Nsec;
        imuData(j).Time(1) = double(nSecs)/1e9;
    end
    
    for i = 2:nMsgs
        msgArray = imuMsgs(i).Mimus;
        
        nIMUs = numel(msgArray);
        
        for j = 1:nIMUs
            msg = msgArray(j).Imu;
            
            acc = [msg.LinearAcceleration.X, msg.LinearAcceleration.Y, msg.LinearAcceleration.Z];
            omega = [msg.AngularVelocity.X, msg.AngularVelocity.Y, msg.AngularVelocity.Z];
            
            % Fill time
            secs = msgArray(j).Imu.Header.Stamp.Sec;
            nSecs = msgArray(j).Imu.Header.Stamp.Nsec;
            t_ij = double(secs - sec0) + double(nSecs)/1e9;
            
            if ismember(t_ij,imuData(j).Time)
                % Missing message
                imuData(j).Time(i) = imuData(j).Time(i-1) + dtIMU;
                imuData(j).RawData(i,:) = nan(1,6);
            else
                % Fill ith element
                imuData(j).Time(i) = t_ij;
                imuData(j).RawData(i,:) = [acc, omega];
            end
        end
    end

    dataset.IMUData = imuData;
end

%% Read mocap ROS data
topics = rosBagSelect.AvailableTopics.Row;
mocapTopics = topics(contains(topics,'vrpn_client_node'));

mocapData = MotiveDataset();

nTrackers = numel(mocapTopics);
trackers(1:nTrackers) = MotiveRigidBody();

for i = 1:nTrackers
    topic = mocapTopics{i};
    mocapSelect = rosBagSelect.select('Topic',topic);
    mocapMsgs = readMessages(mocapSelect,'DataFormat','struct');
    mocapMsgs = [mocapMsgs{:}];
    
    trackerName = regexp(topic,'/vrpn_client_node/(?<trackerName>.+)/pose','names');
    trackerName = trackerName.trackerName;
    
    nMsgs = numel(mocapMsgs);
    
    S = struct();
    
    if isempty(sec0)
        sec0 = mocapMsgs(1).Header.Stamp.Sec;
    end
    
    tj = zeros(nMsgs,1);
    pj = zeros(nMsgs,3);
    rj = zeros(nMsgs,4);
    
    for j = 1:nMsgs
        msg = mocapMsgs(j);
        
        secs = msg.Header.Stamp.Sec;
        nSecs = msg.Header.Stamp.Nsec;
        
        tj(j) = double(secs - sec0) + double(nSecs)/1e9;
        
        pos = [msg.Pose.Position.X, msg.Pose.Position.Y, msg.Pose.Position.Z];
        q = [msg.Pose.Orientation.X, msg.Pose.Orientation.Y,...
            msg.Pose.Orientation.Z, msg.Pose.Orientation.W];
        
        pj(j,:) = pos;
        rj(j,:) = q;
    end
    
    trackers(i).Time = tj;
    trackers(i).Position = pj;
    trackers(i).Rotation = rj;
    trackers(i).Name = trackerName;
end

mocapData.ROSData = trackers;

dataset.MocapData = mocapData;

%% Read robot data
robotData = RobotDataset();

if ismember('/iiwa/joint_states',availableTopics)
    robotSelect = rosBagSelect.select('Topic','/iiwa/joint_states');
    robotMsgs = readMessages(robotSelect,'DataFormat','struct');
    robotMsgs = [robotMsgs{:}];
    
    if isempty(sec0)
        sec0 = robotMsgs(1).Header.Stamp.Sec;
    end
    
    headers = [robotMsgs.Header];
    stamps = [headers.Stamp];
    
    secs = [stamps.Sec].';
    nSecs = [stamps.Nsec].';
    
    robotData.Time = double(secs - sec0) + double(nSecs)/1e9;
    robotData.Position = [robotMsgs.Position].';
    robotData.Velocity = [robotMsgs.Velocity].';
    robotData.Effort = [robotMsgs.Effort].';
    
    dataset.RobotData = robotData;
end

%% Read force data

forceData = ForceDataset();

if ismember('/force',availableTopics)
    forceSelect = rosBagSelect.select('Topic','/force');
    forceMsgs = readMessages(forceSelect,'DataFormat','struct');
    forceMsgs = [forceMsgs{:}];
    
    tForce = forceSelect.MessageList.Time;
    
    if isempty(sec0)
        sec0 = floor(tForce(1));
    end
    
    tForce = tForce - double(sec0);
    
    forceData.Time = tForce;
    forceData.Data = [forceMsgs.Data].';
    
    dataset.ForceData = forceData;
end

end