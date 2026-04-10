function [IEKF, vSensor] = IEKF_init(q_init, u_init, tau_init, Covariances, P0)
% IEKF_INIT Initialize the Iterated Extended Kalman Filter struct and virtual sensor struct.
%   IMU-only version: measurement vector y = [angvel; linacc] for each IMU (36x1 total).
%   State vector x = [q; u; tau] where q = joint angles, u = joint velocities, tau = joint torques.
%
%   Inputs:
%       q_init      - (ND x 1) Initial joint angles [rad]
%       u_init      - (ND x 1) Initial joint velocities [rad/s]
%       tau_init    - (ND x 1) Initial joint torques [Nm]
%       Covariances - Struct with fields Noise_Cov_Gyro_k, Noise_Cov_Acc_k (k=1..6)
%       P0          - (Dx x Dx) Initial state covariance matrix
%
%   Outputs:
%       IEKF    - IEKF filter struct with all storage arrays initialized
%       vSensor - Virtual IMU sensor struct for storing predicted IMU outputs

    ND = length(q_init);
    x_init = [q_init; u_init; tau_init];
    Dx = numel(x_init);
    P_init = P0;

    assert(all(size(x_init) == [Dx, 1]));
    assert(all(size(P_init) == [Dx, Dx]));

    %% IEKF Struct
    IEKF = struct;
    IEKF.Dx       = Dx;
    IEKF.q_init   = q_init;
    IEKF.u_init   = u_init;
    IEKF.tau_init = tau_init;

    % State history arrays (initialized empty, grown during filtering)
    IEKF.x_preds    = NaN(Dx, 0);
    IEKF.P_preds    = NaN(Dx, Dx, 0);
    IEKF.x_upds     = NaN(Dx, 0);
    IEKF.P_upds     = NaN(Dx, Dx, 0);
    IEKF.x_updsIter = NaN(Dx, 0);

    % Jacobian history
    IEKF.F_prev = NaN(Dx, Dx, 0);
    IEKF.L_prev = NaN(Dx, Dx, 0);

    % Torque random-walk standard deviations
    stdTorque = 0.5;
    IEKF.stdJointTorque = stdTorque * ones(ND, 1);

    % Process noise covariance Q: zero noise on q and u, random walk on tau
    q_noise   = zeros(1, ND);
    u_noise   = zeros(1, ND);
    tau_noise = (stdTorque^2) * ones(1, ND);
    IEKF.Q = diag([q_noise, u_noise, tau_noise]);

    % Measurement noise covariance R: block diagonal of scaled IMU covariances
    gyroScales = [150, 250, 250, 250, 250, 150];
    accScales  = [50,  50,  100, 50,  150, 150];
    imuCovBlocks = cell(1, 2 * ND);
    for k = 1:ND
        imuCovBlocks{2*k-1} = Covariances.(['Noise_Cov_Gyro_' num2str(k)]) * gyroScales(k);
        imuCovBlocks{2*k}   = Covariances.(['Noise_Cov_Acc_'  num2str(k)]) * accScales(k);
    end
    IEKF.R = blkdiag(imuCovBlocks{:});

    % Measurement dimension: 6 IMUs x 6 outputs (3 angvel + 3 linacc)
    Dy = size(IEKF.R, 1);
    IEKF.Dy = Dy;

    % Measurement update history arrays
    IEKF.H_current    = NaN(Dy, Dx, 0);
    IEKF.M_current    = NaN(Dy, Dy, 0);
    IEKF.S_current    = NaN(Dy, Dy, 0);
    IEKF.K_current    = NaN(Dx, Dy, 0);
    IEKF.e_current    = NaN(Dy, 0);
    IEKF.h_current    = NaN(Dy, 0);
    IEKF.hat_y_localF = NaN(Dy, 0);

    IEKF.T = [];

    % Initialize at t=0 with zeros (no measurement yet)
    IEKF.H_current(:,:,end+1) = zeros(Dy, Dx);
    IEKF.M_current(:,:,end+1) = zeros(Dy, Dy);
    IEKF.S_current(:,:,end+1) = zeros(Dy, Dy);
    IEKF.K_current(:,:,end+1) = zeros(Dx, Dy);
    IEKF.e_current(:,end+1)    = zeros(Dy, 1);
    IEKF.h_current(:,end+1)    = zeros(Dy, 1);
    IEKF.hat_y_localF(:,end+1) = zeros(Dy, 1);

    % Set initial prediction and update to initialized state at t=0
    IEKF.x_preds(:,end+1)    = x_init;
    IEKF.P_preds(:,:,end+1)  = P_init;
    IEKF.x_upds(:,end+1)     = x_init;
    IEKF.P_upds(:,:,end+1)   = P_init;
    IEKF.x_updsIter(:,end+1) = x_init;
    IEKF.T(end+1) = 0;

    %% Virtual Sensor Struct (stores predicted IMU outputs at each time step)
    vSensor = struct;
    for k = 1:ND
        ks = num2str(k);
        vSensor.(['IMU' ks '_AngVel'])  = zeros(3, 1);
        vSensor.(['IMU' ks '_LinAcc'])  = zeros(3, 1);
        vSensor.(['IMU' ks '_Gravity']) = zeros(3, 1);
    end
    vSensor.T = 0;

end
