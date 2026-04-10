classdef MotiveDataset < AbstractObject
    % MOTIVEDATASET - Class that contains motion capture data from Motive and
    % handles the synchronization between the Motive data and the data from a
    % ROS VRPN node.
    %
    properties (Description = 'Data')
        MotiveData MotiveTracker
        ROSData MotiveTracker
    end
    
    properties (Dependent, Description = 'Data')
        Trackers
        Time
    end
    
    properties (Description = 'Input')
        MotiveFile
    end
    
    %% Get methods
    methods
        function trackerNames = get.Trackers(obj)
            namesA = {obj.MotiveData.Name};
            namesB = {obj.ROSData.Name};
            trackerNames = union(namesA,namesB);
        end
        
        function t = get.Time(obj)
            if ~isempty(obj.MotiveData)
                t = obj.MotiveData(1).Time;
            else
                t = [];
            end
        end
    end
    
    %% Public methods
    methods
        function obj = MotiveDataset(varargin)
            % MotiveDataset constructor
            obj = obj@AbstractObject(varargin{:});
        end
        
        function parse(obj,varargin)
            obj.MotiveData = readMotiveExport(obj.MotiveFile,varargin{:});
        end
        
        function syncMotiveExport(obj,T1,T2)
            % Find a RigidBody tracker to sync
            % trackerNames = setdiff(obj.MotiveData.Properties.VariableNames,'Time');
            allTrackers = obj.MotiveData;
            rigidBodyIdx = arrayfun(@(t) isa(t,'MotiveRigidBody'),allTrackers);
            rigidBodies = allTrackers(rigidBodyIdx);
            
            
            
            for i = 1:numel(rigidBodies)
                rosTracker = obj.ROSData.select(rigidBodies(i).Name);
                if ~isempty(rosTracker)
                    rosTracker = rosTracker(1);
                    rosPos = rosTracker.Position;
                    rawPos = rigidBodies(i).Position;
                    rosTime = rosTracker.Time;
                    rawTime = rigidBodies(i).Time;
                    break
                end
            end
            
            idx_sync_raw = rawTime > T1(1) & rawTime < T1(2);
            idx_sync_ros = rosTime > T2(1) & rosTime < T2(2);
            
            rosPos_sync = rosPos(idx_sync_ros,:);
            rosTime_sync = rosTime(idx_sync_ros);
            rawPos_sync = rawPos(idx_sync_raw,:);
            rawTime_sync = rawTime(idx_sync_raw);
            
            % Take as the correlation signal the distance from the starting point
            rosDist = sqrt(sum((rosPos_sync - rosPos_sync(1,:)).^2,2));
            rawDist = sqrt(sum((rawPos_sync - rawPos_sync(1,:)).^2,2));
            
            fig = figure();
            subplot(2,2,1);
            plot(rawTime_sync,rawDist);
            title('Motive data');
            
            subplot(2,2,3);
            plot(rosTime_sync,rosDist);
            title('ROS data');
            
            [corr,lags] = xcorr(rosDist,rawDist);
            [maxCorr,idxMaxCorr] = max(corr);
            lagMaxCorr = lags(idxMaxCorr);
            
            dt = mean(diff(rawTime));
            
            Nsync = min(numel(rawTime_sync),numel(rosTime_sync));
            ii = (1-lagMaxCorr):(Nsync);
            ii = ii(ii>0);
            ii = ii(ii < Nsync - lags(idxMaxCorr));
            delay = mean(rawTime_sync(ii) - rosTime_sync(ii+lags(idxMaxCorr)));
            
            subplot(2,2,2);
            plot(lags*dt,corr); hold on;
            xline(lagMaxCorr*dt,'k-','LineWidth',2);
            text(lagMaxCorr*dt,maxCorr/2,sprintf('Delay = %.2fs',delay),'BackgroundColor','w')
            title('Correlation');
            
            ButtonH = uicontrol('Style','pushbutton','String','OK','Units','normalized','Position',[0.65 0.2 0.1 0.1],'Visible','on');
            ButtonH.Callback = @removeDelay;
            
            function removeDelay(~,~,~)
                for j = 1:numel(obj.MotiveData)
                    obj.MotiveData(j).Time = obj.MotiveData(j).Time - delay;
                end
                fprintf('Synchronized.\n');
                close(fig);
            end
        end
        
        function Y = getData(obj,idx)
            if nargin > 1 && ~isempty(idx)
                if iscell(idx) || ischar(idx) || isstring(idx)
                    names = {obj.MotiveData.Name};
                    [~,idx] = ismember(idx,names);
                end
            else
                idx = 1:numel(obj.MotiveData);
            end

            Y = [obj.MotiveData(idx).Position];
            
        end
    end
    
end