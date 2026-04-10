function [IEKF, vSensor] = IEKF_update_step_(IEKF, vSensor, Xsens, IMU2, IMU3, IMU4, ...
    osimModel, state, t, torque2Gen, torque3Gen, torque4Gen, ground, GravityVec, ...
    Q1Set, Q2Set, Q3Set, ii, TData, Zero_T_update)
% IEKF_UPDATE_STEP_ Performs one IEKF measurement update iteration.
%
%   Computes virtual IMU measurements from OpenSim, assembles the measurement
%   Jacobian H via numerical perturbation, calculates the Kalman gain, and
%   updates the state estimate.
%
%   Inputs:
%       IEKF       - IEKF struct
%       vSensor    - Virtual sensor struct
%       Xsens      - Struct with real IMU measurements
%       IMU2, IMU3, IMU4 - OpenSim PhysicalOffsetFrame handles for IMUs
%       osimModel  - OpenSim Model
%       state      - OpenSim State
%       t          - Current time step index
%       torque2Gen, torque3Gen, torque4Gen - Torque actuator handles
%       ground     - OpenSim Ground frame
%       GravityVec - Gravity vector as Vec3
%       Q1Set, Q2Set, Q3Set - Coordinate handles
%       ii         - IEKF iteration number
%       TData      - Zero-torque measurement data struct
%       Zero_T_update - Time index after which zero-torque updates are applied

    import org.opensim.modeling.*;

    % Retrieve previous state estimates
    x_prev = IEKF.x_upds(:, end);
    x_iter = IEKF.x_updsIter(:, end);
    P_prev = IEKF.P_upds(:,:, end);

    % Set OpenSim state to current iteration values
    setOpenSimState(Q1Set, Q2Set, Q3Set, torque2Gen, torque3Gen, torque4Gen, state, x_iter);
    osimModel.realizeAcceleration(state);

    % Store virtual torque sensor readings
    vSensor.T1(:, end+1) = torque2Gen.getValue;
    vSensor.T2(:, end+1) = torque3Gen.getValue;
    vSensor.T3(:, end+1) = torque4Gen.getValue;

    % --- Compute predicted IMU measurements at the operating point ---
    [angvel2, linacc2, grav2] = getIMUMeasurements(IMU2, ground, state, GravityVec);
    [angvel3, linacc3, grav3] = getIMUMeasurements(IMU3, ground, state, GravityVec);
    [angvel4, linacc4, grav4] = getIMUMeasurements(IMU4, ground, state, GravityVec);

    % Store in vSensor
    vSensor.IMU2_AngVel(:, end+1)  = angvel2;
    vSensor.IMU2_LinAcc(:, end+1)  = linacc2;
    vSensor.IMU2_Gravity(:, end+1) = grav2;
    vSensor.IMU3_AngVel(:, end+1)  = angvel3;
    vSensor.IMU3_LinAcc(:, end+1)  = linacc3;
    vSensor.IMU3_Gravity(:, end+1) = grav3;
    vSensor.IMU4_AngVel(:, end+1)  = angvel4;
    vSensor.IMU4_LinAcc(:, end+1)  = linacc4;
    vSensor.IMU4_Gravity(:, end+1) = grav4;

    % Predicted measurement vector h(x_iter): [gyro2; acc2+g; gyro3; acc3+g; gyro4; acc4+g; tau]
    h_op = [angvel2; linacc2 + grav2; ...
            angvel3; linacc3 + grav3; ...
            angvel4; linacc4 + grav4; ...
            vSensor.T1(:, end); vSensor.T2(:, end); vSensor.T3(:, end)];

    Ny = length(h_op);
    Nx = IEKF.Dx;

    % --- Numerical Jacobian H via forward finite differences ---
    perturbStep = 1e-8;
    statePerturbed = State(state);
    J = zeros(Ny, Nx);

    for i = 1:Nx
        OPp = x_iter;
        OPp(i) = OPp(i) + perturbStep;

        setOpenSimState(Q1Set, Q2Set, Q3Set, torque2Gen, torque3Gen, torque4Gen, statePerturbed, OPp);
        osimModel.realizeAcceleration(statePerturbed);

        % Perturbed IMU measurements
        [av2p, la2p, g2p] = getIMUMeasurements(IMU2, ground, statePerturbed, GravityVec);
        [av3p, la3p, g3p] = getIMUMeasurements(IMU3, ground, statePerturbed, GravityVec);
        [av4p, la4p, g4p] = getIMUMeasurements(IMU4, ground, statePerturbed, GravityVec);

        h_perturbed = [av2p; la2p + g2p; ...
                       av3p; la3p + g3p; ...
                       av4p; la4p + g4p; ...
                       torque2Gen.getValue; torque3Gen.getValue; torque4Gen.getValue];

        J(:, i) = (h_perturbed - h_op) / perturbStep;
    end

    J = round(J, 7);
    H = J;
    IEKF.H_current(:,:, end+1) = H;

    % Measurement noise Jacobian M = I (noise added directly)
    IEKF.M_current(:,:, end+1) = eye(Ny);

    % Innovation covariance S = H * P * H' + R
    IEKF.S_current(:,:, end+1) = H * P_prev * H' + IEKF.R;

    % Kalman gain K = P * H' / S
    IEKF.K_current(:,:, end+1) = P_prev * H' / IEKF.S_current(:,:, end);

    % Store predicted measurement
    IEKF.h_current(:,:, end+1) = h_op;

    % Linearized predicted measurement: h(x_iter) + H * (x_prev - x_iter)
    IEKF.hat_y_localF(:,:, end+1) = h_op + H * (x_prev - x_iter);

    % --- Assemble actual measurement vector ---
    y_meas = [Xsens.IMU2_AngVel(:, t+1); Xsens.IMU2_LinAcc(:, t+1); ...
              Xsens.IMU3_AngVel(:, t+1); Xsens.IMU3_LinAcc(:, t+1); ...
              Xsens.IMU4_AngVel(:, t+1); Xsens.IMU4_LinAcc(:, t+1); ...
              TData.T1(:, t+1);  TData.T2(:, t+1);  TData.T3(:, t+1)];

    % Apply zero-torque update: tighten torque measurement noise after release
    if t > Zero_T_update
        IEKF.R(19,19) = 0.04;
        IEKF.R(20,20) = 0.01;
        IEKF.R(21,21) = 0.01;
    end

    % Innovation (measurement residual)
    IEKF.e_current(:,:, end+1) = y_meas - IEKF.hat_y_localF(:, end);

    % State update: x_upd = x_prev + K * e
    x_update = x_prev + IEKF.K_current(:,:, end) * IEKF.e_current(:,:, end);
    IEKF.x_updsIter(:, end+1) = x_update;
end

% =========================================================================
% Helper: Set OpenSim state from a 9x1 state vector [q; u; tau]
% =========================================================================
function setOpenSimState(Q1Set, Q2Set, Q3Set, torque2Gen, torque3Gen, torque4Gen, state, x)
    Q1Set.setValue(state, x(1));
    Q2Set.setValue(state, x(2));
    Q3Set.setValue(state, x(3));
    Q1Set.setSpeedValue(state, x(4));
    Q2Set.setSpeedValue(state, x(5));
    Q3Set.setSpeedValue(state, x(6));
    torque2Gen.setValue(x(7));
    torque3Gen.setValue(x(8));
    torque4Gen.setValue(x(9));
end

% =========================================================================
% Helper: Extract angular velocity, linear acceleration, and gravity
%         from an OpenSim IMU frame, all expressed in the local sensor frame.
% =========================================================================
function [angvel, linacc, gravity] = getIMUMeasurements(imuFrame, ground, state, GravityVec)
    % Angular velocity in sensor frame
    gAngVel = imuFrame.getAngularVelocityInGround(state);
    lAngVel = ground.expressVectorInAnotherFrame(state, gAngVel, imuFrame);
    angvel  = [lAngVel.get(0); lAngVel.get(1); lAngVel.get(2)];

    % Linear acceleration in sensor frame
    gLinAcc = imuFrame.getLinearAccelerationInGround(state);
    lLinAcc = ground.expressVectorInAnotherFrame(state, gLinAcc, imuFrame);
    linacc  = [lLinAcc.get(0); lLinAcc.get(1); lLinAcc.get(2)];

    % Gravity vector in sensor frame
    lGrav   = ground.expressVectorInAnotherFrame(state, GravityVec, imuFrame);
    gravity = [lGrav.get(0); lGrav.get(1); lGrav.get(2)];
end
