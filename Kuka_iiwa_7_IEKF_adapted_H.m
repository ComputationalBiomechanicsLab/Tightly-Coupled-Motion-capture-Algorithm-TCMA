%% Kuka_iiwa_7_IEKF_adapted_H.m
% Tightly coupled IMU-based motion capture for the KUKA LBR iiwa 7 R800.
% Estimates joint kinematics and kinetics via an Iterated Extended Kalman Filter (IEKF)
% coupled with an OpenSim multibody dynamic model.

clear; close all; clc;

%% Configuration
SampleFreq = 100;                   % IMU sample rate [Hz]
startSample = 50;                   % ~0.5 s into the recording
endSample   = 2000;                 % ~20 s into the recording
geometryPath = 'C:\OpenSim 4.4\Geometry';
modelFilename = 'Kukaiiwa7_Markers_IMUs.osim';
NUM_IMUS = 6;
NUM_MARKERS_PER_BODY = 4;

%% Import OpenSim libraries and load model
import org.opensim.modeling.*;
ModelVisualizer.addDirToGeometrySearchPaths(geometryPath);
osimModel = Model(modelFilename);

%% Load data files
load MARKERSLOCATIONIMUMARKERS_USE4.mat;
load Rfilter2;        R1 = R;
load 20211103-fast@100Hz_transformed.mat;
load processQ;
load Pnode;
load IMU_markerDatakuka10;
%% Extract IMU data from estimationData struct
imuAccData  = cell(1, NUM_IMUS);
imuGyroData = cell(1, NUM_IMUS);
for k = 1:NUM_IMUS
    imuAccData{k}  = estimationData.IMUData(1, k).LinearAcceleration(:,1:3)';
    imuGyroData{k} = estimationData.IMUData(1, k).AngularVelocity(:,1:3)';
end

%% Interpolate NaNs in IMU data
allIMU = [imuAccData, imuGyroData];
for i = 1:length(allIMU)
    arr = allIMU{i};
    for row = 1:size(arr, 1)
        nanIdx = isnan(arr(row, :));
        if ~any(nanIdx), continue; end
        goodIdx = find(~nanIdx);
        if length(goodIdx) >= 2
            arr(row, nanIdx) = interp1(goodIdx, arr(row, goodIdx), find(nanIdx), 'spline', 'extrap');
        elseif length(goodIdx) == 1
            arr(row, nanIdx) = arr(row, goodIdx);
        end
    end
    allIMU{i} = arr;
end
imuAccData  = allIMU(1:NUM_IMUS);
imuGyroData = allIMU(NUM_IMUS+1:end);

%% Map IMU indices to body segments (IMU4->body1, IMU5->body2, ..., IMU1->body6)
bodyToIMU = [4, 5, 6, 3, 2, 1];
Xsens = struct;
for body = 1:NUM_IMUS
    src = bodyToIMU(body);
    Xsens.(['IMU' num2str(body) '_LinAcc']) = imuAccData{src};
    Xsens.(['IMU' num2str(body) '_AngVel']) = imuGyroData{src};
end

%% Extract MoCap marker data
MocaP = struct;
for k = 1:NUM_IMUS * NUM_MARKERS_PER_BODY
    body = ceil(k / NUM_MARKERS_PER_BODY);
    marker = mod(k - 1, NUM_MARKERS_PER_BODY) + 1;
    fname = sprintf('IMU_Frame%d_Marker%d', body, marker);
    MocaP.(fname) = estimationData.MocapData.MotiveData(1, k).Position';
end

%% Prepend zero column and trim to desired time window
imuFields = fieldnames(Xsens);
for i = 1:length(imuFields)
    Xsens.(imuFields{i}) = [zeros(3,1), Xsens.(imuFields{i})];
end
markerFields = fieldnames(MocaP);
for i = 1:length(markerFields)
    MocaP.(markerFields{i}) = [zeros(3,1), MocaP.(markerFields{i})];
end

%% Create time vector
dTime = 1 / SampleFreq;
nSamples = size(Xsens.IMU1_AngVel, 2) - 2;
Xsens.T = (0:nSamples) * dTime;
Xsens.TotalSimulationTimeValue = length(Xsens.T);

%% Initial joint angles from robot encoders
q_init = estimationData.RobotData.Position(1, 1:NUM_IMUS)';

%% Enable visualizer
osimModel.setUseVisualizer(true);

%% Create coordinate actuators and controllers for all joints
coordSet = osimModel.getCoordinateSet();
optimalForce = 1;
ND = NUM_IMUS;

rotCoords   = cell(1, ND);
actuators   = cell(1, ND);
brains      = cell(1, ND);
torqueGens  = cell(1, ND);

for j = 1:ND
    rotCoords{j} = coordSet.get(j-1);

    actuators{j} = CoordinateActuator(char(rotCoords{j}.getName));
    actuators{j}.setOptimalForce(optimalForce);
    actuators{j}.setName(rotCoords{j}.getName());
    actuators{j}.setMaxControl(Inf);
    actuators{j}.setMinControl(-Inf);
    osimModel.addForce(actuators{j});

    brains{j} = PrescribedController();
    brains{j}.setName(['Brain' num2str(j)]);
    brains{j}.addActuator(actuators{j});

    rotCoords{j}.setDefaultValue(q_init(j));
    rotCoords{j}.setDefaultSpeedValue(0);

    torqueGens{j} = Constant(0);
    brains{j}.prescribeControlForActuator(string(actuators{j}), torqueGens{j});
    osimModel.addController(brains{j});
end

%% Collect initial state values
q_init = zeros(ND, 1);
u_init = zeros(ND, 1);
for j = 1:ND
    q_init(j) = rotCoords{j}.getDefaultValue;
    u_init(j) = rotCoords{j}.getDefaultSpeedValue;
end
load tau_init.mat;
tau_init = tau_init(1:ND, :);

%% Set gravity
ground = osimModel.getGround();
Gravity = -9.81;
osimModel.setGravity(Vec3(0, Gravity, 0));
GravityVec = Vec3(0, -Gravity, 0);

%% Register markers and IMU frames on the model
linkNames = arrayfun(@(k) ['link_' num2str(k)], 1:ND, 'UniformOutput', false);
Links     = cell(1, ND);
IMUFrames = cell(1, ND);
MarkerHandles = cell(ND, NUM_MARKERS_PER_BODY);

for k = 1:ND
    Links{k} = osimModel.getBodySet.get(linkNames{k});

    % Register 4 markers per body
    for m = 1:NUM_MARKERS_PER_BODY
        mName = sprintf('IMUF%dM%d', k, m);
        xField = sprintf('IMU_Frame%dMarker%dX', k, m);
        yField = sprintf('IMU_Frame%dMarker%dY', k, m);
        zField = sprintf('IMU_Frame%dMarker%dZ', k, m);

        mkr = Marker();
        mkr.setName(mName);
        mkr.setParentFrame(Links{k});
        mkr.set_location(Vec3(Orientation.(xField), Orientation.(yField), Orientation.(zField)));
        osimModel.addMarker(mkr);
        MarkerHandles{k, m} = mkr;
    end

    % Register IMU frame
    imuFrame = PhysicalOffsetFrame();
    imuNames = {'FirstIMU','SecondIMU','ThirdIMU','FourthIMU','FifthIMU','SixthIMU'};
    imuFrame.setName(imuNames{k});
    imuFrame.setParentFrame(Links{k});
    txField = sprintf('Xsens_%d_Trans_X', k);
    tyField = sprintf('Xsens_%d_Trans_Y', k);
    tzField = sprintf('Xsens_%d_Trans_Z', k);
    rxField = sprintf('Xsens_%d_Rot_X', k);
    ryField = sprintf('Xsens_%d_Rot_Y', k);
    rzField = sprintf('Xsens_%d_Rot_Z', k);
    imuFrame.set_translation(Vec3(Orientation.(txField), Orientation.(tyField), Orientation.(tzField)));
    imuFrame.set_orientation(Vec3(Orientation.(rxField), Orientation.(ryField), Orientation.(rzField)));
    Links{k}.addComponent(imuFrame);
    IMUFrames{k} = imuFrame;
end

load('Covariances6IMUs.mat');
load('covmarker3_3.mat');

%% Initialize IEKF
[IEKF, vSensor] = IEKF_init(q_init, u_init, tau_init, Covariances, Qaug, cov_marker, R1, P0);

%% Set up reporter for IMU outputs
reportTimeInterval = dTime;
reporter = TableReporterVec3();
reporter.setName('reporter');
reporter.set_report_time_interval(reportTimeInterval);
for k = 1:ND
    reporter.addToReport(IMUFrames{k}.getOutput('linear_acceleration'), ['IMU' num2str(k) '_lin_acc']);
    reporter.addToReport(IMUFrames{k}.getOutput('angular_velocity'),    ['IMU' num2str(k) '_ang_vel']);
end

body7 = osimModel.getBodySet().get('link_7');
osimModel.addComponent(reporter);

%% Initialize system
osimModel.finalizeConnections();
state = osimModel.initSystem();

sviz = osimModel.updVisualizer().updSimbodyVisualizer();
sviz.setShowSimTime(true);

%% Save configured model
osimModel.print('Kukaiiwa7_Markers_IMUs_configured.osim');
fprintf('Model printed.\n');

%% Run IEKF estimation loop
tic

IEKF.dTime = reportTimeInterval;
IEKF.finalTime = Xsens.T(end);
IEKF.n = round(IEKF.finalTime / IEKF.dTime);
IEKF.NoI = 3;

% Cache coordinate set handle for performance
cachedCoordSet = osimModel.getCoordinateSet();

h_wait = waitbar(0, 'Running IEKF...');
manager = Manager(osimModel);
manager.initialize(state);

coords = zeros(IEKF.n + 1, 3);

for t = 1:IEKF.n

    % --- IEKF Predict Step Part 1: set prior states ---
    IEKF = IEKF_predict_step_part_1(IEKF, osimModel, state, ...
        torqueGens{1}, torqueGens{2}, torqueGens{3}, ...
        torqueGens{4}, torqueGens{5}, torqueGens{6}, t);

    % Restore torques after Jacobian perturbations
    for j = 1:ND
        torqueGens{j}.setValue(IEKF.tau_prev(j, 1));
    end
    osimModel.realizeAcceleration(state);

    % --- IEKF Predict Step Manager: integrate EoM ---
    state.setTime((t-1) * IEKF.dTime);
    state = manager.integrate(t * IEKF.dTime);

    % --- IEKF Predict Step Part 2: get predicted states ---
    IEKF = IEKF_predict_step_part_2(IEKF, osimModel, state, t);

    % --- IEKF Update Step ---

    IEKF.H_current1 = 0;
    IEKF.M_current1 = 0;
    IEKF.K_current1 = 0;
    max_iter = 3;
    conv_tol = 1e-6;
    
     for ii = 1:max_iter
         if t>1
        x_iter_prev = x_iter;
         else 
            x_iter_prev = 0; 
         end
    [x_iter, IEKF, vSensor] = IEKF_update_step_adapted_H(IEKF, MocaP, vSensor, Xsens, ...
        IMUFrames{1}, IMUFrames{2}, IMUFrames{3}, IMUFrames{4}, IMUFrames{5}, IMUFrames{6}, ...
        MarkerHandles{1,1}, MarkerHandles{1,2}, MarkerHandles{1,3}, MarkerHandles{1,4}, ...
        MarkerHandles{2,1}, MarkerHandles{2,2}, MarkerHandles{2,3}, MarkerHandles{2,4}, ...
        MarkerHandles{3,1}, MarkerHandles{3,2}, MarkerHandles{3,3}, MarkerHandles{3,4}, ...
        MarkerHandles{4,1}, MarkerHandles{4,2}, MarkerHandles{4,3}, MarkerHandles{4,4}, ...
        MarkerHandles{5,1}, MarkerHandles{5,2}, MarkerHandles{5,3}, MarkerHandles{5,4}, ...
        MarkerHandles{6,1}, MarkerHandles{6,2}, MarkerHandles{6,3}, MarkerHandles{6,4}, ...
        osimModel, state, t, ...
        torqueGens{1}, torqueGens{2}, torqueGens{3}, ...
        torqueGens{4}, torqueGens{5}, torqueGens{6}, ...
        ground, GravityVec);
       dx = x_iter - x_iter_prev;
        if norm(dx) / max(norm(x_iter), 1) < conv_tol
            break;
        end
     end
    % --- IEKF Step 3: finalize update ---
    IEKF.x_upds(:, end) = IEKF.x_updsIter(:, end);

    P_prev = IEKF.P_upds(:, :, end);
    P_update = (eye(IEKF.Dx) - IEKF.K_current1 * IEKF.H_current1) * P_prev;
    IEKF.P_upds(:, :, end) = P_update;

    % Set model to updated states
    for j = 1:ND
        cachedCoordSet.get(j-1).setValue(state, IEKF.x_upds(j, end));
        cachedCoordSet.get(j-1).setSpeedValue(state, IEKF.x_upds(ND + j, end));
        torqueGens{j}.setValue(IEKF.x_upds(2*ND + j, end));
    end
    osimModel.realizeAcceleration(state);

    % Track end-effector position
    transform7 = body7.getTransformInGround(state);
    pos7 = transform7.p();
    coords(t+1, :) = [pos7.get(0), pos7.get(1), pos7.get(2)];

    waitbar(t / IEKF.n, h_wait);
end
close(h_wait);

elapsed = toc;
fprintf('IEKF completed in %.1f seconds.\n', elapsed);
