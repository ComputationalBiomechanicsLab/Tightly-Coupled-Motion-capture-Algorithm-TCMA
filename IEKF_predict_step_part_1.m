function IEKF = IEKF_predict_step_part_1(IEKF, osimModel, state, torque1Gen, torque2Gen, torque3Gen, torque4Gen, torque5Gen, torque6Gen, t)
% IEKF_PREDICT_STEP_PART_1 Set prior states and compute Jacobians F_prev and L_prev.
%   Computes the continuous-time Jacobian G = dg/dx via numerical perturbation,
%   then discretizes via matrix exponential: F = expm(G * dt).

    import org.opensim.modeling.*;

    ND = length(IEKF.q_init);
    Dx = 3 * ND; % [q; u; tau]

    % Extract previous updated states
    IEKF.q_prev   = IEKF.x_upds(1:ND, end);
    IEKF.u_prev   = IEKF.x_upds(ND+1:2*ND, end);
    IEKF.tau_prev = IEKF.x_upds(2*ND+1:end, end);

    % Set OpenSim model to previous state
    coordSet = osimModel.getCoordinateSet();
    for j = 1:ND
        coordSet.get(j-1).setValue(state, IEKF.q_prev(j));
        coordSet.get(j-1).setSpeedValue(state, IEKF.u_prev(j));
    end
    torqueGens = {torque1Gen, torque2Gen, torque3Gen, torque4Gen, torque5Gen, torque6Gen};
    for j = 1:ND
        torqueGens{j}.setValue(IEKF.tau_prev(j));
    end

    % Realize accelerations to get QDot and UDot at the operating point
    osimModel.realizeAcceleration(state);

    % Collect operating point values
    OP = [IEKF.q_prev; IEKF.u_prev; IEKF.tau_prev];

    qDotStates = state.getQDot;
    uDotStates = state.getUDot;
    QDotOP = zeros(ND, 1);
    UDotOP = zeros(ND, 1);
    for j = 1:ND
        QDotOP(j) = qDotStates.get(j-1);
        UDotOP(j) = uDotStates.get(j-1);
    end
    TauDotOP = zeros(ND, 1); % tau_dot = 0 (random walk model)
    EvalOP = [QDotOP; UDotOP; TauDotOP];

    % Numerical Jacobian via perturbation
    h_pert = 1e-8;
    OPperturb = OP;
    J = NaN(Dx, Dx);

    for i = 1:Dx
        OPperturb(i) = OPperturb(i) + h_pert;

        statePerturbed = State(state);
        for j = 1:ND
            coordSet.get(j-1).setValue(statePerturbed, OPperturb(j));
            coordSet.get(j-1).setSpeedValue(statePerturbed, OPperturb(ND + j));
            torqueGens{j}.setValue(OPperturb(2*ND + j));
        end

        osimModel.realizeAcceleration(statePerturbed);

        qDotP = statePerturbed.getQDot;
        uDotP = statePerturbed.getUDot;
        for j = 1:ND
            J(j, i)      = (round(qDotP.get(j-1), 15) - round(QDotOP(j), 15)) / h_pert;
            J(ND+j, i)   = (round(uDotP.get(j-1), 15) - round(UDotOP(j), 15)) / h_pert;
            J(2*ND+j, i) = 0; % tau_dot = 0, no dependency on states
        end

        OPperturb(i) = OP(i); % Reset perturbation
    end

    J = round(J, 7);
    F_discrete = expm(J * IEKF.dTime);
    IEKF.F_prev(:, :, end+1) = F_discrete;

    % L_prev: identity (additive process noise assumption)
    L_discrete = expm(eye(Dx) * IEKF.dTime);
    IEKF.L_prev(:, :, end+1) = L_discrete;

    % Restore torques to operating point values
    for j = 1:ND
        torqueGens{j}.setValue(IEKF.tau_prev(j));
    end

end
