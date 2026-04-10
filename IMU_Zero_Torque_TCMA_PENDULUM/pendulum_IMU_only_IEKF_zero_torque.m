% PENDULUM_IMU_ONLY_IEKF_ZERO_TORQUE
% Main script for TCMA_{I,0T}: IMU-only IEKF with zero-torque updates
% on a 3-DoF pendulum system.
%
% This script:
%   1. Loads the OpenSim model and IMU sensor data
%   2. Registers IMU frames and markers on each link
%   3. Creates coordinate actuators with prescribed controllers
%   4. Initializes the IEKF
%   5. Runs the predict-update loop over all time steps
%   6. Saves results

clear all; close all; clc;
import org.opensim.modeling.*;

%% ======================== Model Setup ========================
path = 'C:\OpenSim 4.4\Geometry';
ModelVisualizer.addDirToGeometrySearchPaths(path);

osimModel = Model('Pendulum.osim');
osimModel.setUseVisualizer(true);

%% ======================== Load IMU Data ========================
SampleFreq = 128;  % Hz
dTime = 1 / SampleFreq;

IMUz = readtable('imuData01.xls');
IMU11 = table2array(IMUz);
IMU11(1, :) = 0;  % Zero the first sample

Xsens = struct;
Xsens.IMU2_AngVel = IMU11(:, 20:22)';  % Upper link
Xsens.IMU2_LinAcc = IMU11(:, 17:19)';
Xsens.IMU3_AngVel = IMU11(:, 35:37)';  % Middle link
Xsens.IMU3_LinAcc = IMU11(:, 32:34)';
Xsens.IMU4_AngVel = IMU11(:, 50:52)';  % Lower link
Xsens.IMU4_LinAcc = IMU11(:, 47:49)';

nSamples = size(IMU11, 1);

% Zero-torque data (torques are zero throughout for this passive system)
TData = struct;
TData.T1 = zeros(1, nSamples);
TData.T2 = zeros(1, nSamples);
TData.T3 = zeros(1, nSamples);

load('cov.mat');

%% ======================== Select Time Window ========================
StartTime = 0.1;
FinalTime = 8;

startIdx = StartTime * SampleFreq;
endIdx   = FinalTime * SampleFreq;

% Trim and prepend a zero column for t=0
Xsens.IMU2_AngVel = [zeros(3,1), Xsens.IMU2_AngVel(:, startIdx:endIdx)];
Xsens.IMU2_LinAcc = [zeros(3,1), Xsens.IMU2_LinAcc(:, startIdx:endIdx)];
Xsens.IMU3_AngVel = [zeros(3,1), Xsens.IMU3_AngVel(:, startIdx:endIdx)];
Xsens.IMU3_LinAcc = [zeros(3,1), Xsens.IMU3_LinAcc(:, startIdx:endIdx)];
Xsens.IMU4_AngVel = [zeros(3,1), Xsens.IMU4_AngVel(:, startIdx:endIdx)];
Xsens.IMU4_LinAcc = [zeros(3,1), Xsens.IMU4_LinAcc(:, startIdx:endIdx)];

TData.T1 = [0, TData.T1(:, startIdx:endIdx)];
TData.T2 = [0, TData.T2(:, startIdx:endIdx)];
TData.T3 = [0, TData.T3(:, startIdx:endIdx)];

% Detect pendulum release: first zero-crossing after the dominant peak
TH1 = max(Xsens.IMU4_LinAcc(1, :)) * 0.75;
[~, peak_locs] = findpeaks(Xsens.IMU4_LinAcc(1, :), 'MinPeakHeight', TH1);
y = Xsens.IMU4_LinAcc(1, peak_locs(1):end);
zero_crossings = find(diff(sign(y)));
Zero_T_update = zero_crossings(1) + peak_locs(1);

% Build time vector
nSteps = size(Xsens.IMU2_AngVel, 2) - 1;  % Exclude the prepended zero
Xsens.T = (0:nSteps-1) * dTime;
Xsens.TotalSimulationTimeValue = length(Xsens.T);

%% ======================== Coordinate Actuators ========================
coordSet = osimModel.getCoordinateSet();
optimalForce = 1;

% Joint 1 (upper)
Rotation1_2 = coordSet.get(6);
Rotation1_2Actuator = CoordinateActuator(char(Rotation1_2.getName));
Rotation1_2Actuator.setOptimalForce(optimalForce);
Rotation1_2Actuator.setName(Rotation1_2.getName());
Rotation1_2Actuator.setMaxControl(Inf);
Rotation1_2Actuator.setMinControl(-Inf);
osimModel.addForce(Rotation1_2Actuator);

% Joint 2 (middle)
Rotation2_3 = coordSet.get(7);
Rotation2_3Actuator = CoordinateActuator(char(Rotation2_3.getName));
Rotation2_3Actuator.setOptimalForce(optimalForce);
Rotation2_3Actuator.setName(Rotation2_3.getName());
Rotation2_3Actuator.setMaxControl(Inf);
Rotation2_3Actuator.setMinControl(-Inf);
osimModel.addForce(Rotation2_3Actuator);

% Joint 3 (lower)
Rotation3_4 = coordSet.get(8);
Rotation3_4Actuator = CoordinateActuator(char(Rotation3_4.getName));
Rotation3_4Actuator.setOptimalForce(optimalForce);
Rotation3_4Actuator.setName(Rotation3_4.getName());
Rotation3_4Actuator.setMaxControl(Inf);
Rotation3_4Actuator.setMinControl(-Inf);
osimModel.addForce(Rotation3_4Actuator);

%% ======================== Gravity ========================
ground = osimModel.getGround();
Gravity = -9.81;
osimModel.setGravity(Vec3(0, Gravity, 0));
GravityVec = Vec3(0, -Gravity, 0);

%% ======================== Initial Conditions ========================
q_init = [0; 0; 0.08618436];

% Set default values on coordinates
Rotation1_2.setDefaultValue(q_init(1));
Rotation1_2.setDefaultSpeedValue(0);
Rotation2_3.setDefaultValue(q_init(2));
Rotation2_3.setDefaultSpeedValue(0);
Rotation3_4.setDefaultValue(q_init(3));
Rotation3_4.setDefaultSpeedValue(0);

% Prescribed controllers with initial zero torques
act2 = string(Rotation1_2Actuator);
act3 = string(Rotation2_3Actuator);
act4 = string(Rotation3_4Actuator);

brain2 = PrescribedController();
brain2.setName('Brain2');
brain2.addActuator(Rotation1_2Actuator);
torque2Gen = Constant(0);
brain2.prescribeControlForActuator(act2, torque2Gen);
osimModel.addController(brain2);

brain3 = PrescribedController();
brain3.setName('Brain3');
brain3.addActuator(Rotation2_3Actuator);
torque3Gen = Constant(0);
brain3.prescribeControlForActuator(act3, torque3Gen);
osimModel.addController(brain3);

brain4 = PrescribedController();
brain4.setName('Brain4');
brain4.addActuator(Rotation3_4Actuator);
torque4Gen = Constant(0);
brain4.prescribeControlForActuator(act4, torque4Gen);
osimModel.addController(brain4);

u_init   = [0; 0; 0];
tau_init = [0; 0; 0];

%% ======================== Register IMUs & Markers ========================
load('Orient9.mat');

% --- Static link markers ---
Link1 = osimModel.getBodySet.get('static');
addMarkerToModel(osimModel, 'StaticO', Link1, Orientation, 'StaticO');
addMarkerToModel(osimModel, 'StaticX', Link1, Orientation, 'StaticX');
addMarkerToModel(osimModel, 'StaticY', Link1, Orientation, 'StaticY');

% --- Upper link: IMU2 + markers ---
Link2 = osimModel.getBodySet.get('upper');
IMU2 = PhysicalOffsetFrame();
IMU2.setName('upper_imu');
IMU2.setParentFrame(Link2);
IMU2.set_translation(Vec3(Orientation.Xsens_2_Trans_X, Orientation.Xsens_2_Trans_Y, Orientation.Xsens_2_Trans_Z));
IMU2.set_orientation(Vec3(Orientation.Xsens_2_Rot_X, Orientation.Xsens_2_Rot_Y, Orientation.Xsens_2_Rot_Z));
Link2.addComponent(IMU2);

addMarkerToModel(osimModel, 'UpperO',  Link2, Orientation, 'UpperO');
addMarkerToModel(osimModel, 'UpperX',  Link2, Orientation, 'UpperX');
addMarkerToModel(osimModel, 'UpperY',  Link2, Orientation, 'UpperY');
addMarkerToModel(osimModel, 'UpperJ1', Link2, Orientation, 'UpperJ1');

% --- Middle link: IMU3 + markers ---
Link3 = osimModel.getBodySet.get('middle');
IMU3 = PhysicalOffsetFrame();
IMU3.setName('middle_imu');
IMU3.setParentFrame(Link3);
IMU3.set_translation(Vec3(Orientation.Xsens_3_Trans_X, Orientation.Xsens_3_Trans_Y, Orientation.Xsens_3_Trans_Z));
IMU3.set_orientation(Vec3(Orientation.Xsens_3_Rot_X, Orientation.Xsens_3_Rot_Y, Orientation.Xsens_3_Rot_Z));
Link3.addComponent(IMU3);

addMarkerToModel(osimModel, 'MiddleO',  Link3, Orientation, 'MiddleO');
addMarkerToModel(osimModel, 'MiddleX',  Link3, Orientation, 'MiddleX');
addMarkerToModel(osimModel, 'MiddleY',  Link3, Orientation, 'MiddleY');
addMarkerToModel(osimModel, 'MiddleJ1', Link3, Orientation, 'MiddleJ1');

% --- Lower link: IMU4 + markers ---
Link4 = osimModel.getBodySet.get('lower');
IMU4 = PhysicalOffsetFrame();
IMU4.setName('lower_imu');
IMU4.setParentFrame(Link4);
IMU4.set_translation(Vec3(Orientation.Xsens_4_Trans_X, Orientation.Xsens_4_Trans_Y, Orientation.Xsens_4_Trans_Z));
IMU4.set_orientation(Vec3(Orientation.Xsens_4_Rot_X, Orientation.Xsens_4_Rot_Y, Orientation.Xsens_4_Rot_Z));
Link4.addComponent(IMU4);

addMarkerToModel(osimModel, 'LowerO', Link4, Orientation, 'LowerO');
addMarkerToModel(osimModel, 'LowerX', Link4, Orientation, 'LowerX');
addMarkerToModel(osimModel, 'LowerY', Link4, Orientation, 'LowerY');
addMarkerToModel(osimModel, 'Lower5', Link4, Orientation, 'Lower5');
addMarkerToModel(osimModel, 'Lower6', Link4, Orientation, 'Lower6');

%% ======================== Initialize IEKF ========================
[IEKF, vSensor] = IEKF_init(q_init, u_init, tau_init, Covariances);

reportTimeInterval = dTime;
reporter = TableReporterVec3();
reporter.setName('reporter');
reporter.set_report_time_interval(reportTimeInterval);
osimModel.addComponent(reporter);

%% ======================== Initialize OpenSim System ========================
state = osimModel.initSystem();
sviz = osimModel.updVisualizer().updSimbodyVisualizer();
sviz.setShowSimTime(true);
osimModel.print('Pendulum12.osim');

%% ======================== IEKF Main Loop ========================
tic

IEKF.dTime     = reportTimeInterval;
IEKF.finalTime = Xsens.T(end);
IEKF.n         = round(IEKF.finalTime / IEKF.dTime);
IEKF.NoI       = 3;  % Max IEKF iterations (used 2 in the loop below)

h_waitbar = waitbar(0, 'Running TCMA...');
manager   = Manager(osimModel);
manager.initialize(state);

Q1Set = osimModel.getCoordinateSet().get(6);
Q2Set = osimModel.getCoordinateSet().get(7);
Q3Set = osimModel.getCoordinateSet().get(8);

% Preallocate energy arrays
KinEnergy    = zeros(1, IEKF.n);
PotEnergy    = zeros(1, IEKF.n);

for t = 1:IEKF.n
    % --- Prediction step ---
    IEKF = IEKF_predict_step_part_1(IEKF, osimModel, state, ...
        torque2Gen, torque3Gen, torque4Gen, t, Q1Set, Q2Set, Q3Set);

    % Restore torques after Jacobian perturbation
    torque2Gen.setValue(IEKF.tau_prev(1));
    torque3Gen.setValue(IEKF.tau_prev(2));
    torque4Gen.setValue(IEKF.tau_prev(3));
    osimModel.realizeAcceleration(state);

    % Integrate equations of motion over one time step
    state.setTime((t-1) * IEKF.dTime);
    state = manager.integrate(t * IEKF.dTime);

    IEKF = IEKF_predict_step_part_2(IEKF, osimModel, state, t);

    % --- Measurement update (2 IEKF iterations) ---
    for ii = 1:2
        [IEKF, vSensor] = IEKF_update_step_(IEKF, vSensor, Xsens, ...
            IMU2, IMU3, IMU4, osimModel, state, t, ...
            torque2Gen, torque3Gen, torque4Gen, ground, GravityVec, ...
            Q1Set, Q2Set, Q3Set, ii, TData, Zero_T_update);
    end

    % Finalize update: set iterated estimate as the update
    IEKF.x_upds(:, end) = IEKF.x_updsIter(:, end);

    % Joseph-form covariance update: P = (I - K*H) * P_pred
    P_prev   = IEKF.P_upds(:,:, end);
    P_update = (eye(IEKF.Dx) - IEKF.K_current(:,:, end) * IEKF.H_current(:,:, end)) * P_prev;
    IEKF.P_upds(:,:, end) = P_update;

    % Apply updated state to OpenSim
    Q1Set.setValue(state, IEKF.x_upds(1, end));
    Q2Set.setValue(state, IEKF.x_upds(2, end));
    Q3Set.setValue(state, IEKF.x_upds(3, end));
    Q1Set.setSpeedValue(state, IEKF.x_upds(4, end));
    Q2Set.setSpeedValue(state, IEKF.x_upds(5, end));
    Q3Set.setSpeedValue(state, IEKF.x_upds(6, end));
    torque2Gen.setValue(IEKF.x_upds(7, end));
    torque3Gen.setValue(IEKF.x_upds(8, end));
    torque4Gen.setValue(IEKF.x_upds(9, end));
    osimModel.realizeAcceleration(state);

    % Record energy
    KinEnergy(t) = osimModel.calcKineticEnergy(state);
    PotEnergy(t) = osimModel.calcPotentialEnergy(state);

    waitbar(t / IEKF.n, h_waitbar);
end
close(h_waitbar);

elapsedTime = toc;
fprintf('TCMA completed in %.1f seconds (%.1fx real-time)\n', ...
    elapsedTime, elapsedTime / IEKF.finalTime);

%% ======================== Save Results ========================
IEKF_IMU_ZeroTorque = IEKF;

% =========================================================================
% Local function: Add a marker to the OpenSim model
% =========================================================================
function addMarkerToModel(osimModel, markerName, parentFrame, Orientation, prefix)
    import org.opensim.modeling.*;
    m = Marker();
    m.setName(markerName);
    m.setParentFrame(parentFrame);
    tx = Orientation.([prefix, '_Trans_X']);
    ty = Orientation.([prefix, '_Trans_Y']);
    tz = Orientation.([prefix, '_Trans_Z']);
    m.set_location(Vec3(tx, ty, tz));
    osimModel.addMarker(m);
end
