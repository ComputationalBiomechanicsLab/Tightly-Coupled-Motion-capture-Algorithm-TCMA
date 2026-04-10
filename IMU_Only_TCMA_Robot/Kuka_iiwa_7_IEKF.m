%% Kuka_iiwa_7_IEKF.m
% IMU-only motion estimation for the KUKA LBR iiwa 7 R800.
% Estimates joint kinematics and kinetics via an Iterated Extended Kalman Filter (IEKF)
% coupled with an OpenSim multibody dynamic model. Uses only IMU data (no markers).

clear; close all; clc;

%% Configuration
SampleFreq    = 100;                              % IMU sample rate [Hz]
geometryPath  = 'C:\OpenSim 4.4\Geometry';
modelFilename = 'Kukaiiwa7_IMUs.osim';
NUM_IMUS      = 6;

%% Import OpenSim libraries and load model
import org.opensim.modeling.*;
ModelVisualizer.addDirToGeometrySearchPaths(geometryPath);
osimModel = Model(modelFilename);

%% Load data files
load MARKERSLOCATIONIMUMARKERS_USE4.mat;
load processQ;
load IMU_markerDatakuka10;
load 20211103-fast@100Hz_transformed.mat;

%% Extract IMU data from estimationData struct
imuAccData  = cell(1, NUM_IMUS);
imuGyroData = cell(1, NUM_IMUS);
for k = 1:NUM_IMUS
    imuAccData{k}  = estimationData.IMUData(1, k).CalibratedData(:, 1:3)';
    imuGyroData{k} = estimationData.IMUData(1, k).CalibratedData(:, 4:6)';
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

%% Prepend zero column for t=0 placeholder
imuFields = fieldnames(Xsens);
for i = 1:length(imuFields)
    Xsens.(imuFields{i}) = [zeros(3, 1), Xsens.(imuFields{i})];
end

%% Create time vector for the IMU measured signals
dTime = 1 / SampleFreq;
nSamples = size(Xsens.IMU1_AngVel, 2) - 1;
Xsens.T = (0:nSamples-1) * dTime;
Xsens.TotalSimulationTimeValue = length(Xsens.T);

%% Visualization prompt
fprintf('\n');
str1 = input('Visualize the simulation? [y]es/[n]o: ', 's');
if str1 == 'y'
    osimModel.setUseVisualizer(true);
    fprintf('Visualization enabled.\n');
end

%% Get ground reference and set gravity
ground = osimModel.getGround();
Gravity = -9.81;
osimModel.setGravity(Vec3(0, Gravity, 0));
GravityVec = Vec3(0, -Gravity, 0); % Accelerometer measures reaction force

%% Set initial conditions from robot encoder data
q_init_raw = estimationData.RobotData.Position(1, 1:6)';
PerturbValue = deg2rad(0);

%% Create coordinate actuators and controllers for all joints
coordSet = osimModel.getCoordinateSet();
ND = NUM_IMUS;
optimalForce = 1;

torqueGens = cell(1, ND);
q_angles   = zeros(ND, 1);
for k = 1:ND
    coord = coordSet.get(k-1);

    % Create actuator
    actuator = CoordinateActuator(char(coord.getName));
    actuator.setOptimalForce(optimalForce);
    actuator.setName(coord.getName());
    actuator.setMaxControl(Inf);
    actuator.setMinControl(-Inf);
    osimModel.addForce(actuator);

    % Set initial angle (with optional perturbation for joints 2-6)
    if k == 1
        q_angles(k) = q_init_raw(k);
    else
        q_angles(k) = q_init_raw(k) + PerturbValue;
    end
    coord.setDefaultValue(q_angles(k));
    coord.setDefaultSpeedValue(deg2rad(0));

    % Create controller with constant torque
    brain = PrescribedController();
    brain.setName(['Brain' num2str(k)]);
    brain.addActuator(actuator);
    torqueGens{k} = Constant(0);
    brain.prescribeControlForActuator(string(actuator), torqueGens{k});
    osimModel.addController(brain);
end

% Collect initial state vectors
q_init = zeros(ND, 1);
u_init = zeros(ND, 1);
for k = 1:ND
    q_init(k) = coordSet.get(k-1).getDefaultValue;
    u_init(k) = coordSet.get(k-1).getDefaultSpeedValue;
end
load tau_init.mat;
tau_init = tau_init(1:ND, :);

%% Create virtual IMU frames attached to robot links
Links     = cell(1, ND);
IMUFrames = cell(1, ND);
imuNames  = {'FirstIMU', 'SecondIMU', 'ThirdIMU', 'FourthIMU', 'FifthIMU', 'SixthIMU'};

for k = 1:ND
    Links{k} = osimModel.getBodySet.get(['link_' num2str(k)]);

    imuFrame = PhysicalOffsetFrame();
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

%% Load covariances and initialize IEKF
load('Covariances6IMUs.mat');
load Pnode;
[IEKF, vSensor] = IEKF_init(q_init, u_init, tau_init, Covariances, P0);

%% Set up reporter for IMU outputs
reportTimeInterval = dTime;
reporter = TableReporterVec3();
reporter.setName('reporter');
reporter.set_report_time_interval(reportTimeInterval);
for k = 1:ND
    reporter.addToReport(IMUFrames{k}.getOutput('linear_acceleration'), ['IMU' num2str(k) '_lin_acc']);
    reporter.addToReport(IMUFrames{k}.getOutput('angular_velocity'),    ['IMU' num2str(k) '_ang_vel']);
end
osimModel.addComponent(reporter);

%% Initialize system
osimModel.finalizeConnections();
state = osimModel.initSystem();

if str1 == 'y'
    sviz = osimModel.updVisualizer().updSimbodyVisualizer();
    sviz.setShowSimTime(true);
end

%% Save configured model
osimModel.print('Kukaiiwa7_IMUs_configured.osim');
fprintf('Model printed.\n');

%% Run IEKF estimation loop
tic

IEKF.dTime     = reportTimeInterval;
IEKF.finalTime = Xsens.T(end);
IEKF.n         = round(IEKF.finalTime / IEKF.dTime);
IEKF.NoI       = 3;

% Cache coordinate set handle for performance
cachedCoordSet = osimModel.getCoordinateSet();

h_wait = waitbar(0, 'Running IEKF...');
manager = Manager(osimModel);
manager.initialize(state);

for t = 1:IEKF.n

    % --- IEKF Predict Step Part 1: set prior states and compute Jacobians ---
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

    % --- IEKF Update Step (iterated) ---
    IEKF.H_current1 = 0;
    IEKF.M_current1 = 0;
    IEKF.K_current1 = 0;
    for i = 1:2
        [IEKF, vSensor] = IEKF_update_step(IEKF, vSensor, Xsens, ...
            IMUFrames{1}, IMUFrames{2}, IMUFrames{3}, ...
            IMUFrames{4}, IMUFrames{5}, IMUFrames{6}, ...
            osimModel, state, t, ...
            torqueGens{1}, torqueGens{2}, torqueGens{3}, ...
            torqueGens{4}, torqueGens{5}, torqueGens{6}, ...
            ground, GravityVec);
    end

    % --- IEKF Step 3: finalize update ---
    IEKF.x_upds(:, end) = IEKF.x_updsIter(:, end);

    P_prev   = IEKF.P_upds(:, :, end);
    P_update = (eye(IEKF.Dx) - IEKF.K_current1 * IEKF.H_current1) * P_prev;
    IEKF.P_upds(:, :, end) = P_update;

    % Set model to updated states
    for j = 1:ND
        cachedCoordSet.get(j-1).setValue(state, IEKF.x_upds(j, end));
        cachedCoordSet.get(j-1).setSpeedValue(state, IEKF.x_upds(ND + j, end));
        torqueGens{j}.setValue(IEKF.x_upds(2*ND + j, end));
    end
    osimModel.realizeAcceleration(state);

    waitbar(t / IEKF.n, h_wait);
end
close(h_wait);

elapsed = toc;
fprintf('IEKF completed in %.1f seconds.\n', elapsed);
