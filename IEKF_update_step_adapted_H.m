function [x_iter, IEKF, vSensor] = IEKF_update_step_adapted_H(IEKF, MocaP, vSensor, Xsens, ...
    IMU1, IMU2, IMU3, IMU4, IMU5, IMU6, ...
    IMUF1M1, IMUF1M2, IMUF1M3, IMUF1M4, ...
    IMUF2M1, IMUF2M2, IMUF2M3, IMUF2M4, ...
    IMUF3M1, IMUF3M2, IMUF3M3, IMUF3M4, ...
    IMUF4M1, IMUF4M2, IMUF4M3, IMUF4M4, ...
    IMUF5M1, IMUF5M2, IMUF5M3, IMUF5M4, ...
    IMUF6M1, IMUF6M2, IMUF6M3, IMUF6M4, ...
    osimModel, state, t, ...
    torque1Gen, torque2Gen, torque3Gen, torque4Gen, torque5Gen, torque6Gen, ...
    ground, GravityVec)
% IEKF_UPDATE_STEP_ADAPTED_H Measurement update step of the IEKF.
%   Computes H Jacobian via numerical perturbation, Kalman gain, and state update.

    import org.opensim.modeling.*;

    ND = length(IEKF.q_init);
    x_prev = IEKF.x_upds(:, end);
    x_iter = IEKF.x_updsIter(:, end);
    P_prev = IEKF.P_upds(:, :, end);

    % Set model state to x_iter
    coordSet = osimModel.getCoordinateSet();
    for j = 1:ND
        coordSet.get(j-1).setValue(state, x_iter(j));
        coordSet.get(j-1).setSpeedValue(state, x_iter(ND+j));
    end
    torqueGens = {torque1Gen, torque2Gen, torque3Gen, torque4Gen, torque5Gen, torque6Gen};
    for j = 1:ND
        torqueGens{j}.setValue(x_iter(2*ND+j));
    end

    osimModel.realizeAcceleration(state);

    %% Evaluate h(x_iter): virtual sensor outputs at operating point
    IMUs = {IMU1, IMU2, IMU3, IMU4, IMU5, IMU6};
    markerHandles = {IMUF1M1, IMUF1M2, IMUF1M3, IMUF1M4, ...
                     IMUF2M1, IMUF2M2, IMUF2M3, IMUF2M4, ...
                     IMUF3M1, IMUF3M2, IMUF3M3, IMUF3M4, ...
                     IMUF4M1, IMUF4M2, IMUF4M3, IMUF4M4, ...
                     IMUF5M1, IMUF5M2, IMUF5M3, IMUF5M4, ...
                     IMUF6M1, IMUF6M2, IMUF6M3, IMUF6M4};

    % Get virtual IMU readings and marker positions at operating point
    imuAngVelOP = zeros(3, ND);
    imuLinAccOP = zeros(3, ND); % linacc + gravity (what the accelerometer measures)
    markerPosOP = zeros(3, 24);

    for k = 1:ND
        % Angular velocity in local IMU frame
        gAngVel = IMUs{k}.getAngularVelocityInGround(state);
        lAngVel = ground.expressVectorInAnotherFrame(state, gAngVel, IMUs{k});
        localAV = [lAngVel.get(0); lAngVel.get(1); lAngVel.get(2)];
        imuAngVelOP(:, k) = localAV;
        vSensor.(['IMU' num2str(k) '_AngVel'])(:, end+1) = localAV;

        % Gravity in local IMU frame
        lGrav = ground.expressVectorInAnotherFrame(state, GravityVec, IMUs{k});
        localG = [lGrav.get(0); lGrav.get(1); lGrav.get(2)];
        vSensor.(['IMU' num2str(k) '_Gravity'])(:, end+1) = localG;

        % Linear acceleration in local IMU frame (acc + gravity = what accelerometer reads)
        gLinAcc = IMUs{k}.getLinearAccelerationInGround(state);
        lLinAcc = ground.expressVectorInAnotherFrame(state, gLinAcc, IMUs{k});
        localLA = [lLinAcc.get(0); lLinAcc.get(1); lLinAcc.get(2)];
        imuLinAccOP(:, k) = localLA + localG; % Accelerometer measures linacc + gravity
        vSensor.(['IMU' num2str(k) '_LinAcc'])(:, end+1) = localLA;
    end

    % Marker positions in ground frame
    for m = 1:24
        loc = markerHandles{m}.getLocationInGround(state);
        markerPosOP(:, m) = [loc.get(0); loc.get(1); loc.get(2)];
    end
    % Store marker positions in vSensor
    markerNames = {'IMUF1M1','IMUF1M2','IMUF1M3','IMUF1M4', ...
                   'IMUF2M1','IMUF2M2','IMUF2M3','IMUF2M4', ...
                   'IMUF3M1','IMUF3M2','IMUF3M3','IMUF3M4', ...
                   'IMUF4M1','IMUF4M2','IMUF4M3','IMUF4M4', ...
                   'IMUF5M1','IMUF5M2','IMUF5M3','IMUF5M4', ...
                   'IMUF6M1','IMUF6M2','IMUF6M3','IMUF6M4'};
    for m = 1:24
        vSensor.(markerNames{m}) = markerPosOP(:, m);
    end

    % Build EvaluateOP: [angvel1; linacc1; angvel2; linacc2; ...; markerPos1; markerPos2; ...]
    EvaluateOP = zeros(36 + 72, 1); % 6 IMUs * 6 + 24 markers * 3
    for k = 1:ND
        idx = (k-1)*6;
        EvaluateOP(idx+1:idx+3) = imuAngVelOP(:, k);
        EvaluateOP(idx+4:idx+6) = imuLinAccOP(:, k);
    end
    for m = 1:24
        idx = 36 + (m-1)*3;
        EvaluateOP(idx+1:idx+3) = markerPosOP(:, m);
    end

    %% Compute measurement Jacobian H via numerical perturbation
    numMeas = length(EvaluateOP);
    numStates = 3 * ND;
    h_pert = 1e-8;
    OPperturb = x_iter;
    J = NaN(numMeas, numStates);

    for i = 1:numStates
        OPperturb(i) = OPperturb(i) + h_pert;

        stateP = State(state);
        for j = 1:ND
            coordSet.get(j-1).setValue(stateP, OPperturb(j));
            coordSet.get(j-1).setSpeedValue(stateP, OPperturb(ND+j));
            torqueGens{j}.setValue(OPperturb(2*ND+j));
        end
        osimModel.realizeAcceleration(stateP);

        % Perturbed IMU readings
        EvalPerturb = zeros(numMeas, 1);
        for k = 1:ND
            gAV = IMUs{k}.getAngularVelocityInGround(stateP);
            lAV = ground.expressVectorInAnotherFrame(stateP, gAV, IMUs{k});

            lGP = ground.expressVectorInAnotherFrame(stateP, GravityVec, IMUs{k});

            gLA = IMUs{k}.getLinearAccelerationInGround(stateP);
            lLA = ground.expressVectorInAnotherFrame(stateP, gLA, IMUs{k});

            idx = (k-1)*6;
            EvalPerturb(idx+1) = lAV.get(0);
            EvalPerturb(idx+2) = lAV.get(1);
            EvalPerturb(idx+3) = lAV.get(2);
            EvalPerturb(idx+4) = lLA.get(0) + lGP.get(0);
            EvalPerturb(idx+5) = lLA.get(1) + lGP.get(1);
            EvalPerturb(idx+6) = lLA.get(2) + lGP.get(2);
        end

        % Perturbed marker positions
        for m = 1:24
            locP = markerHandles{m}.getLocationInGround(stateP);
            idx = 36 + (m-1)*3;
            EvalPerturb(idx+1) = locP.get(0);
            EvalPerturb(idx+2) = locP.get(1);
            EvalPerturb(idx+3) = locP.get(2);
        end

        % Jacobian column
        J(:, i) = (round(EvalPerturb, 15) - round(EvaluateOP, 15)) / h_pert;

        OPperturb(i) = x_iter(i); % Reset
    end
    J = round(J, 7);

    % Restore torques after perturbation
    for j = 1:ND
        torqueGens{j}.setValue(x_iter(2*ND+j));
    end

    %% Assemble actual measurement vector y from sensor data
    % IMU measurements
    y_imu = zeros(36, 1);
    for k = 1:ND
        idx = (k-1)*6;
        y_imu(idx+1:idx+3) = Xsens.(['IMU' num2str(k) '_AngVel'])(:, t+1);
        y_imu(idx+4:idx+6) = Xsens.(['IMU' num2str(k) '_LinAcc'])(:, t+1);
    end

    % Marker measurements (NOTE: IMUF1M4 maps to Frame2_Marker4 per original implementation)
    markerMocapFields = { ...
        'IMU_Frame1_Marker1','IMU_Frame1_Marker2','IMU_Frame1_Marker3','IMU_Frame2_Marker4', ...
        'IMU_Frame2_Marker1','IMU_Frame2_Marker2','IMU_Frame2_Marker3','IMU_Frame2_Marker4', ...
        'IMU_Frame3_Marker1','IMU_Frame3_Marker2','IMU_Frame3_Marker3','IMU_Frame3_Marker4', ...
        'IMU_Frame4_Marker1','IMU_Frame4_Marker2','IMU_Frame4_Marker3','IMU_Frame4_Marker4', ...
        'IMU_Frame5_Marker1','IMU_Frame5_Marker2','IMU_Frame5_Marker3','IMU_Frame5_Marker4', ...
        'IMU_Frame6_Marker1','IMU_Frame6_Marker2','IMU_Frame6_Marker3','IMU_Frame6_Marker4'};
    y_markers = zeros(72, 1);
    for m = 1:24
        y_markers((m-1)*3+1 : m*3) = MocaP.(markerMocapFields{m})(:, t+1);
    end

    y_localX = [y_imu; y_markers];

    %% Handle NaN measurements (missing markers)
    valid_mask = ~isnan(y_localX);
    y_index = find(valid_mask);
    y_localX = y_localX(y_index);
    J = J(y_index, :);

    IEKF.H_current(:, :, end+1) = NaN(size(IEKF.H_current, 1), size(IEKF.H_current, 2)); % placeholder
    H_current = J;
    IEKF.H_current1 = J;

    Ms = length(y_index);
    M_current = eye(Ms);
    IEKF.M_current(:, :, end+1) = NaN(size(IEKF.M_current, 1), size(IEKF.M_current, 2)); % placeholder
    IEKF.M_current1 = M_current;

    %% Compute Kalman gain
    % R submatrix matches the valid measurement dimension (original indexing)
    R2 = IEKF.R(1:Ms, 1:Ms);

    S = H_current * P_prev * H_current' + M_current * R2 * M_current';
    K = P_prev * H_current' / S;
    IEKF.K_current1 = K;

    %% Predicted measurement h(x_iter, 0) and linearized prediction
    h_full = EvaluateOP;
    IEKF.h_current(:, end+1) = h_full;

    h_valid = h_full(y_index);
    hat_y = h_valid + H_current * (x_prev - x_iter);

    %% Innovation and state update
    e = y_localX - hat_y;
    x_update = x_prev + K * e;

    IEKF.x_updsIter(:, end+1) = x_update;

end
