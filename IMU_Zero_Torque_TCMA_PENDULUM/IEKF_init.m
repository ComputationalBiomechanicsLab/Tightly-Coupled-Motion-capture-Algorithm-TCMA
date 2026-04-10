function [IEKF, vSensor] = IEKF_init(q_init, u_init, tau_init, Covariances)
% IEKF_INIT Initializes the IEKF state estimator and virtual sensor structs.
%   [IEKF, vSensor] = IEKF_INIT(q_init, u_init, tau_init, Covariances)
%
%   Inputs:
%       q_init      - Initial generalized coordinates [3x1]
%       u_init      - Initial generalized velocities  [3x1]
%       tau_init    - Initial joint torques            [3x1]
%       Covariances - Struct with sensor noise covariance fields:
%                     Noise_Cov_Gyro_2, Noise_Cov_Acc_2, ..., _3, ..., _4
%
%   Outputs:
%       IEKF    - IEKF struct with state, covariance, and filter matrices
%       vSensor - Virtual sensor struct for predicted IMU measurements

    % Assemble initial state vector x = [q; u; tau]
    x_init = [q_init; u_init; tau_init];
    Dx = numel(x_init);
    P_init = eye(Dx);

    % --- Build IEKF struct ---
    IEKF = struct;
    IEKF.Dx = Dx;
    IEKF.q_init   = q_init;
    IEKF.u_init   = u_init;
    IEKF.tau_init = tau_init;

    % Torque process noise standard deviations
    IEKF.stdJointTorque2 = 1;
    IEKF.stdJointTorque3 = 1;
    IEKF.stdJointTorque4 = 1;

    % Process noise covariance Q: only torque states have nonzero noise
    IEKF.Q = diag([0, 0, 0, ...   % q noise (zero)
                   0, 0, 0, ...   % u noise (zero)
                   IEKF.stdJointTorque2^2, ...
                   IEKF.stdJointTorque3^2, ...
                   IEKF.stdJointTorque4^2]);

    % Measurement noise covariance R: [gyro2, acc2, gyro3, acc3, gyro4, acc4, tau1, tau2, tau3]
    scaleFactor = 150;  % Empirical scaling of sensor covariances
    Cov_torque  = 15.5; % Torque measurement noise variance
    IEKF.R = blkdiag(Covariances.Noise_Cov_Gyro_2 * scaleFactor, ...
                     Covariances.Noise_Cov_Acc_2  * scaleFactor, ...
                     Covariances.Noise_Cov_Gyro_3 * scaleFactor, ...
                     Covariances.Noise_Cov_Acc_3  * scaleFactor, ...
                     Covariances.Noise_Cov_Gyro_4 * scaleFactor, ...
                     Covariances.Noise_Cov_Acc_4  * scaleFactor, ...
                     Cov_torque, Cov_torque, Cov_torque);

    % Measurement dimension
    Dy = size(IEKF.R, 1);
    IEKF.Dy = Dy;

    % --- Preallocate filter history arrays ---
    % Using cell arrays or 3D matrices; initialized with the t=0 values
    IEKF.x_preds    = x_init;
    IEKF.P_preds    = P_init;
    IEKF.x_upds     = x_init;
    IEKF.P_upds     = P_init;
    IEKF.x_updsIter = x_init;
    IEKF.F_prev     = zeros(Dx, Dx, 0);
    IEKF.L_prev     = zeros(Dx, Dx, 0);

    % Measurement model Jacobians and filter matrices at t=0 (no measurement)
    IEKF.H_current  = zeros(Dy, Dx);
    IEKF.M_current  = zeros(Dy, Dy);
    IEKF.S_current  = zeros(Dy, Dy);
    IEKF.K_current  = zeros(Dx, Dy);
    IEKF.e_current  = zeros(Dy, 1);
    IEKF.h_current  = zeros(Dy, 1);
    IEKF.hat_y_localF = zeros(Dy, 1);

    % Time vector
    IEKF.T = 0;

    % --- Virtual sensor struct (predicted IMU measurements) ---
    vSensor = struct;
    zeroCol = zeros(3, 1);

    % IMU 2 (upper link)
    vSensor.IMU2_AngVel   = zeroCol;
    vSensor.IMU2_LinAcc   = zeroCol;
    vSensor.IMU2_Gravity  = zeroCol;

    % IMU 3 (middle link)
    vSensor.IMU3_AngVel   = zeroCol;
    vSensor.IMU3_LinAcc   = zeroCol;
    vSensor.IMU3_Gravity  = zeroCol;

    % IMU 4 (lower link)
    vSensor.IMU4_AngVel   = zeroCol;
    vSensor.IMU4_LinAcc   = zeroCol;
    vSensor.IMU4_Gravity  = zeroCol;

    % Virtual torque sensor readings
    vSensor.T1 = 0;
    vSensor.T2 = 0;
    vSensor.T3 = 0;

    % Time
    vSensor.T = 0;
end
