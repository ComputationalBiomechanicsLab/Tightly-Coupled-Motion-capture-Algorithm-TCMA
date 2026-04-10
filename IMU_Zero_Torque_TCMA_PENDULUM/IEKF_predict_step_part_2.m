function IEKF = IEKF_predict_step_part_2(IEKF, osimModel, state, t)
% IEKF_PREDICT_STEP_PART_2 Extracts the predicted state after numerical
%   integration and computes the predicted state covariance.
%
%   This function is called after OpenSim's integrator has advanced the
%   state by one time step. It reads the integrated q and u from the state,
%   combines them with the (constant) tau, and propagates the covariance.
%
%   Inputs:
%       IEKF      - IEKF struct
%       osimModel - OpenSim Model object
%       state     - OpenSim State object (post-integration)
%       t         - Current time step index

    import org.opensim.modeling.*;

    osimModel.realizeAcceleration(state);

    % Read integrated generalized coordinates (indices 7,8,9 for 3 DoFs)
    q = state.getQ;
    q_pred = [q.get(7); q.get(8); q.get(9)];

    % Read integrated generalized velocities (indices 6,7,8)
    u = state.getU;
    u_pred = [u.get(6); u.get(7); u.get(8)];

    % Torques remain constant during prediction (random walk with zero mean)
    tau_pred = IEKF.tau_prev;

    % Assemble predicted state
    x_pred = [q_pred; u_pred; tau_pred];

    % Predicted covariance: P_pred = F * P_upd * F' + L * Q * L'
    F = IEKF.F_prev(:,:,end);
    L = IEKF.L_prev(:,:,end);
    P_prev = IEKF.P_upds(:,:,end);
    P_pred = F * P_prev * F' + L * IEKF.Q * L';

    % Store predicted state and covariance
    IEKF.x_preds(:, end+1)    = x_pred;
    IEKF.P_preds(:,:, end+1)  = P_pred;
    IEKF.T(end+1) = t * IEKF.dTime;

    % Initialize update arrays for current time step with predicted values
    IEKF.x_upds(:, end+1)     = x_pred;
    IEKF.P_upds(:,:, end+1)   = P_pred;
    IEKF.x_updsIter(:, end+1) = x_pred;
end
