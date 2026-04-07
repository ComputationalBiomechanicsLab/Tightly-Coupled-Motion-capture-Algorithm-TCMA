function IEKF = IEKF_predict_step_part_2(IEKF, osimModel, state, t)
% IEKF_PREDICT_STEP_PART_2 Get predicted states from OpenSim after EoM integration.
%   Computes the predicted state covariance P = F*P*F' + L*Q*L'.

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

    % Initialize update values to prediction
    IEKF.x_upds(:, end+1)     = x;
    IEKF.P_upds(:, :, end+1)  = P;
    IEKF.x_updsIter(:, end+1) = x;

end
