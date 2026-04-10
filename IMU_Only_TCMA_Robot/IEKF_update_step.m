function [IEKF, vSensor] = IEKF_update_step(IEKF, vSensor, Xsens, ...
    IMU1, IMU2, IMU3, IMU4, IMU5, IMU6, ...
    osimModel, state, t, ...
    torque1Gen, torque2Gen, torque3Gen, torque4Gen, torque5Gen, torque6Gen, ...
    ground, GravityVec)
% IEKF_UPDATE_STEP Measurement update step of the IEKF (IMU-only).
%   Computes H Jacobian via numerical perturbation of virtual IMU outputs,
%   then computes innovation covariance S, Kalman gain K, and state update.
%   Measurement vector: y = [angvel_1; linacc_1; angvel_2; linacc_2; ... ; angvel_6; linacc_6] (36x1)
%
%   Inputs:
%       IEKF       - IEKF filter struct
%       vSensor    - Virtual IMU sensor struct
%       Xsens      - Struct with measured IMU data fields IMUk_AngVel, IMUk_LinAcc
%       IMU1..IMU6 - OpenSim PhysicalOffsetFrame handles for virtual IMUs
%       osimModel  - OpenSim Model object
%       state      - OpenSim State object
%       t          - Current time step index
%       torqueXGen - Constant torque generator handles for joints 1-6
%       ground     - OpenSim Ground reference frame
%       GravityVec - Vec3 gravity vector (reaction: [0, +9.81, 0])
%
%   Outputs:
%       IEKF    - Updated struct with H, M, S, K, innovation, and iterated state
%       vSensor - Updated virtual sensor struct with predicted IMU outputs

    import org.opensim.modeling.*;

    ND = length(IEKF.q_init);
    x_prev = IEKF.x_upds(:, end);
    x_iter = IEKF.x_updsIter(:, end);
    P_prev = IEKF.P_upds(:, :, end);

    % Set model state to x_iter
    coordSet = osimModel.getCoordinateSet();
    for j = 1:ND
        coordSet.get(j-1).setValue(state, x_iter(j));
        coordSet.get(j-1).setSpeedValue(state, x_iter(ND + j));
    end
    torqueGens = {torque1Gen, torque2Gen, torque3Gen, torque4Gen, torque5Gen, torque6Gen};
    for j = 1:ND
        torqueGens{j}.setValue(x_iter(2*ND + j));
    end

    osimModel.realizeAcceleration(state);

    %% Evaluate h(x_iter, 0): virtual IMU outputs at operating point
    IMUs = {IMU1, IMU2, IMU3, IMU4, IMU5, IMU6};
    Dy = 6 * ND; % 6 outputs per IMU (3 angvel + 3 linacc)
    EvaluateOP = zeros(Dy, 1);

    for k = 1:ND
        % Angular velocity in local IMU frame
        gAngVel = IMUs{k}.getAngularVelocityInGround(state);
        lAngVel = ground.expressVectorInAnotherFrame(state, gAngVel, IMUs{k});
        localAV = [lAngVel.get(0); lAngVel.get(1); lAngVel.get(2)];
        vSensor.(['IMU' num2str(k) '_AngVel'])(:, end+1) = localAV;

        % Gravity in local IMU frame
        lGrav = ground.expressVectorInAnotherFrame(state, GravityVec, IMUs{k});
        localG = [lGrav.get(0); lGrav.get(1); lGrav.get(2)];
        vSensor.(['IMU' num2str(k) '_Gravity'])(:, end+1) = localG;

        % Linear acceleration in local IMU frame (accelerometer = linacc + gravity)
        gLinAcc = IMUs{k}.getLinearAccelerationInGround(state);
        lLinAcc = ground.expressVectorInAnotherFrame(state, gLinAcc, IMUs{k});
        localLA = [lLinAcc.get(0); lLinAcc.get(1); lLinAcc.get(2)];
        vSensor.(['IMU' num2str(k) '_LinAcc'])(:, end+1) = localLA;

        % Stack into EvaluateOP: [angvel_k; linacc_k + gravity_k]
        idx = (k-1) * 6;
        EvaluateOP(idx+1:idx+3) = localAV;
        EvaluateOP(idx+4:idx+6) = localLA + localG;
    end

    %% Compute measurement Jacobian H via numerical perturbation
    numStates = 3 * ND;
    h_pert = 1e-8;
    OPperturb = x_iter;
    J = NaN(Dy, numStates);

    for i = 1:numStates
        OPperturb(i) = OPperturb(i) + h_pert;

        stateP = State(state);
        for j = 1:ND
            coordSet.get(j-1).setValue(stateP, OPperturb(j));
            coordSet.get(j-1).setSpeedValue(stateP, OPperturb(ND + j));
            torqueGens{j}.setValue(OPperturb(2*ND + j));
        end
        osimModel.realizeAcceleration(stateP);

        % Perturbed IMU readings
        EvalPerturb = zeros(Dy, 1);
        for k = 1:ND
            gAV = IMUs{k}.getAngularVelocityInGround(stateP);
            lAV = ground.expressVectorInAnotherFrame(stateP, gAV, IMUs{k});

            lGP = ground.expressVectorInAnotherFrame(stateP, GravityVec, IMUs{k});

            gLA = IMUs{k}.getLinearAccelerationInGround(stateP);
            lLA = ground.expressVectorInAnotherFrame(stateP, gLA, IMUs{k});

            idx = (k-1) * 6;
            EvalPerturb(idx+1) = lAV.get(0);
            EvalPerturb(idx+2) = lAV.get(1);
            EvalPerturb(idx+3) = lAV.get(2);
            EvalPerturb(idx+4) = lLA.get(0) + lGP.get(0);
            EvalPerturb(idx+5) = lLA.get(1) + lGP.get(1);
            EvalPerturb(idx+6) = lLA.get(2) + lGP.get(2);
        end

        % Jacobian column via finite difference
        J(:, i) = (round(EvalPerturb, 15) - round(EvaluateOP, 15)) / h_pert;

        OPperturb(i) = x_iter(i);
    end

    J = round(J, 7);

    % Restore torques after perturbation
    for j = 1:ND
        torqueGens{j}.setValue(x_iter(2*ND + j));
    end

    % Store H and M Jacobians
    H_current = J;
    IEKF.H_current(:, :, end+1) = H_current;

    M_current = eye(Dy);
    IEKF.M_current(:, :, end+1) = M_current;

    % Temporary storage for main loop covariance update
    IEKF.H_current1 = H_current;
    IEKF.M_current1 = M_current;

    %% Compute innovation covariance and Kalman gain
    S = H_current * P_prev * H_current' + M_current * IEKF.R * M_current';
    IEKF.S_current(:, :, end+1) = S;

    K = P_prev * H_current' / S;
    IEKF.K_current(:, :, end+1) = K;
    IEKF.K_current1 = K;

    %% Store predicted measurement h(x_iter, 0)
    IEKF.h_current(:, end+1) = EvaluateOP;

    %% Compute linearized predicted measurement
    hat_y = EvaluateOP + H_current * (x_prev - x_iter);
    IEKF.hat_y_localF(:, end+1) = hat_y;

    %% Assemble actual measurement vector from IMU data
    y_localX = zeros(Dy, 1);
    for k = 1:ND
        idx = (k-1) * 6;
        y_localX(idx+1:idx+3) = Xsens.(['IMU' num2str(k) '_AngVel'])(:, t+1);
        y_localX(idx+4:idx+6) = Xsens.(['IMU' num2str(k) '_LinAcc'])(:, t+1);
    end

    %% Compute innovation and state update
    e = y_localX - hat_y;
    IEKF.e_current(:, end+1) = e;

    x_update = x_prev + K * e;
    IEKF.x_updsIter(:, end+1) = x_update;

end
