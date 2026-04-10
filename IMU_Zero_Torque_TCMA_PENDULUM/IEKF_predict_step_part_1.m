function IEKF = IEKF_predict_step_part_1(IEKF, osimModel, state, torque2Gen, torque3Gen, torque4Gen, t, Q1Set, Q2Set, Q3Set)
% IEKF_PREDICT_STEP_PART_1 Computes the continuous-time Jacobian F and
%   discretizes it via matrix exponential. Also computes the noise Jacobian L.
%
%   This function:
%     1. Sets the OpenSim state to the previous IEKF-updated values
%     2. Evaluates q-dot and u-dot at the operating point
%     3. Numerically computes the state Jacobian F via forward perturbation
%     4. Discretizes F and L using matrix exponential
%
%   Inputs:
%       IEKF       - IEKF struct (uses x_upds, q_init, u_init, tau_init, dTime)
%       osimModel  - OpenSim Model object
%       state      - OpenSim State object
%       torque2Gen, torque3Gen, torque4Gen - Constant torque actuators
%       t          - Current time step index
%       Q1Set, Q2Set, Q3Set - OpenSim coordinate handles for joints 1-3

    import org.opensim.modeling.*;

    Nq = length(IEKF.q_init);
    Nu = length(IEKF.u_init);
    Ntau = length(IEKF.tau_init);
    Nx = Nq + Nu + Ntau;

    % Extract previous updated states
    IEKF.q_prev   = IEKF.x_upds(1:Nq, end);
    IEKF.u_prev   = IEKF.x_upds(Nq+1:Nq+Nu, end);
    IEKF.tau_prev = IEKF.x_upds(Nq+Nu+1:end, end);

    % Set OpenSim state to previous IEKF estimates
    Q1Set.setValue(state, IEKF.q_prev(1));
    Q2Set.setValue(state, IEKF.q_prev(2));
    Q3Set.setValue(state, IEKF.q_prev(3));
    Q1Set.setSpeedValue(state, IEKF.u_prev(1));
    Q2Set.setSpeedValue(state, IEKF.u_prev(2));
    Q3Set.setSpeedValue(state, IEKF.u_prev(3));
    torque2Gen.setValue(IEKF.tau_prev(1));
    torque3Gen.setValue(IEKF.tau_prev(2));
    torque4Gen.setValue(IEKF.tau_prev(3));

    % Realize accelerations to compute q-dot and u-dot
    osimModel.realizeAcceleration(state);

    % Extract q-dot, u-dot at operating point (indices 6,7,8 are the 3 DoFs)
    qDotStates = state.getQDot;
    uDotStates = state.getUDot;
    EvalOP = [qDotStates.get(6); qDotStates.get(7); qDotStates.get(8); ...
              uDotStates.get(6); uDotStates.get(7); uDotStates.get(8); ...
              0; 0; 0];  % tau-dot = 0 (random walk)

    % Operating point vector
    OP = [IEKF.q_prev; IEKF.u_prev; IEKF.tau_prev];

    % --- Numerical Jacobian via forward finite differences ---
    h = 1e-8;
    statePerturbed = State(state);  % Deep copy for perturbation
    J = zeros(Nx, Nx);

    for i = 1:Nx
        % Perturb the i-th state
        OPp = OP;
        OPp(i) = OPp(i) + h;

        Q1Set.setValue(statePerturbed, OPp(1));
        Q2Set.setValue(statePerturbed, OPp(2));
        Q3Set.setValue(statePerturbed, OPp(3));
        Q1Set.setSpeedValue(statePerturbed, OPp(4));
        Q2Set.setSpeedValue(statePerturbed, OPp(5));
        Q3Set.setSpeedValue(statePerturbed, OPp(6));
        torque2Gen.setValue(OPp(7));
        torque3Gen.setValue(OPp(8));
        torque4Gen.setValue(OPp(9));

        osimModel.realizeAcceleration(statePerturbed);

        qDotP = statePerturbed.getQDot;
        uDotP = statePerturbed.getUDot;
        EvalPerturbed = [qDotP.get(6); qDotP.get(7); qDotP.get(8); ...
                         uDotP.get(6); uDotP.get(7); uDotP.get(8); ...
                         0; 0; 0];

        J(:, i) = (EvalPerturbed - EvalOP) / h;
    end

    J = round(J, 7);

    % Discretize: F_discrete = expm(F_continuous * dt)
    IEKF.F_prev(:,:,end+1) = expm(J * IEKF.dTime);

    % Noise Jacobian L is identity (noise enters all states directly)
    IEKF.L_prev(:,:,end+1) = expm(eye(Nx) * IEKF.dTime);
end
