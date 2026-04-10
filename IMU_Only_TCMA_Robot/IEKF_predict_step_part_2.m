function IEKF = IEKF_predict_step_part_2(IEKF, osimModel, state, t)
% IEKF_PREDICT_STEP_PART_2 Get predicted states from OpenSim after EoM integration.
%   Extracts q, u from the integrated state and computes the predicted state
%   covariance P = F*P_prev*F' + L*Q*L'.
%
%   Inputs:
%       IEKF      - IEKF filter struct
%       osimModel - OpenSim Model object
%       state     - OpenSim State object (after integration)
%       t         - Current time step index
%
%   Output:
%       IEKF      - Updated struct with predicted x, P, and initialized update arrays

    import org.opensim.modeling.*;

    ND = length(IEKF.q_init);

    osimModel.realizeAcceleration(state);

    % Extract predicted states from integrated model
    q_vec = state.getQ;
    u_vec = state.getU;
    q_pred = zeros(ND, 1);
    u_pred = zeros(ND, 1);
    for j = 1:ND
        q_pred(j) = q_vec.get(j-1);
        u_pred(j) = u_vec.get(j-1);
    end

    % Torque prediction: tau_dot = 0, so tau_{k+1} = tau_k
    tau_pred = IEKF.tau_prev;

    x = [q_pred; u_pred; tau_pred];

    % Predicted state covariance: P = F*P_prev*F' + L*Q*L'
    P_prev = IEKF.P_upds(:, :, end);
    F = IEKF.F_prev(:, :, end);
    L = IEKF.L_prev(:, :, end);
    P = F * P_prev * F' + L * IEKF.Q * L';

    % Store predictions
    IEKF.x_preds(:, end+1)    = x;
    IEKF.P_preds(:, :, end+1) = P;
    IEKF.T(end+1)             = t * IEKF.dTime;

    % Initialize update values to prediction (will be overwritten in update step)
    IEKF.x_upds(:, end+1)     = x;
    IEKF.P_upds(:, :, end+1)  = P;
    IEKF.x_updsIter(:, end+1) = x;

end
