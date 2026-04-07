# Tightly Coupled Motion-capture Algorithm (TCMA)

A MATLAB implementation of a tightly coupled multibody dynamics and multi-sensor fusion algorithm for simultaneous kinematics and kinetics estimation, based on an Iterated Extended Kalman Filter (IEKF) coupled with [OpenSim](https://opensim.stanford.edu/) multibody models.

> **Paper:** H. Osman, D. de Kanter, J. Boelens, M. Kok, and A. Seth, "A Tightly Coupled Multibody Dynamics and Multi-Sensor Fusion Algorithm for Simultaneous Kinematics and Kinetics Estimation," 2026.

## Overview

TCMA directly integrates IMU accelerometer and gyroscope measurements with multibody dynamic models to estimate joint kinematics (angles, velocities) and kinetics (torques) simultaneously. Unlike traditional approaches that first compute IMU orientations (requiring magnetometers) and then solve inverse kinematics/dynamics as separate steps, TCMA:

- **Bypasses magnetometer dependence** — uses only accelerometer and gyroscope data, making it suitable for magnetically disturbed environments.
- **Enforces dynamic consistency** — the full equations of motion serve as the process model, ensuring estimated kinematics and kinetics are physically consistent.
- **Supports multi-sensor fusion** — the framework can incorporate optical motion capture markers and virtual zero-torque measurements alongside IMU data.
- **Provides direct kinetic estimates** — eliminates the need for low-pass filtering and numerical differentiation required by conventional inverse dynamics pipelines.

![TCMA Overview](docs/tcma_overview.png)

## Algorithm

The state vector at each time step is:

```
x = [q; q̇; τ]
```

where `q` are joint angles, `q̇` joint velocities, and `τ` joint torques.

**Prediction step:** The state is propagated forward using OpenSim's Runge-Kutta integrator to numerically integrate the multibody equations of motion. Joint torques evolve as a random walk.

**Update step:** An Iterated Extended Kalman Filter refines the predicted state using sensor measurements. The measurement model maps states to expected IMU readings (angular velocity and linear acceleration) and, optionally, marker positions and zero-torque priors. Iteration handles the nonlinearity of the measurement model.

## Repository Structure

```
├── Kuka_iiwa_7_IEKF_adapted_H.m      # Main script: KUKA robot IEKF estimation
├── IEKF_init.m                        # IEKF and virtual sensor initialization
├── IEKF_predict_step_part_1.m         # Prediction step: set prior states & compute Jacobian
├── IEKF_predict_step_part_2.m         # Prediction step: extract predicted states post-integration
├── IEKF_update_step_adapted_H.m       # Update step: iterated measurement update
├── Kukaiiwa7_Markers_IMUs.osim        # OpenSim model (KUKA LBR iiwa 7 R800)
├── Kukaiiwa7_Markers_IMUs_configured.osim  # Configured OpenSim model (with actuators/sensors)
│
├── IMUDataset.m                       # IMU data container class
├── IMUCalibration.m                   # IMU calibration utilities
├── RobotDataset.m                     # Robot encoder/torque data container
├── ForceDataset.m                     # Force data container
├── MotiveDataset.m                    # OptiTrack Motive data container
├── MotiveMarker.m                     # Marker data structure
├── MotiveRigidBody.m                  # Rigid body data structure
├── MotiveTracker.m                    # Tracker data structure
├── ROSDataset.m                       # ROS bag data container
├── readROSBag.m                       # ROS bag reader utility
│
├── DataValidator.m                    # Data validation and synchronization
├── AbstractObject.m                   # Base class for data objects
├── interpolateAndResample.m           # Signal interpolation and resampling
├── cut.m                              # Time-window extraction utility
```

## Requirements

- **MATLAB** R2023b or later (tested on R2025b)
- **OpenSim 4.4+** with the MATLAB scripting interface configured ([installation guide](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab))
- OpenSim Geometry files (set `geometryPath` in the main script)

## Getting Started

1. **Install OpenSim** and configure the MATLAB API following the [official instructions](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab).

2. **Clone this repository:**
   ```bash
   git clone https://github.com/hassanwagih/A-Tightly-Coupled-Multibody-Dynamics-and-Multi-Sensor-Fusion-Algorithm-TCMA-.git
   cd A-Tightly-Coupled-Multibody-Dynamics-and-Multi-Sensor-Fusion-Algorithm-TCMA-
   ```

3. **Prepare your data.** The main script expects `.mat` files containing IMU measurements, marker positions, and robot encoder data in a structured format. See `Kuka_iiwa_7_IEKF_adapted_H.m` for the expected variable names and structure.

4. **Configure paths** in `Kuka_iiwa_7_IEKF_adapted_H.m`:
   ```matlab
   geometryPath  = 'C:\OpenSim 4.4\Geometry';   % Path to OpenSim geometry files
   modelFilename = 'Kukaiiwa7_Markers_IMUs.osim'; % OpenSim model file
   ```

5. **Run the estimation:**
   ```matlab
   Kuka_iiwa_7_IEKF_adapted_H
   ```

## Experimental Validation

TCMA was validated on two physical systems with accurate ground truth:

| System | DoF | Joint Angle RMSD | Joint Torque RMSD |
|--------|-----|-------------------|-------------------|
| 3-DoF Pendulum (IMU only) | 3 | ≤ 3.75° vs. marker IK | ≤ 3.02 Nm vs. marker ID |
| KUKA iiwa Robot (IMU only) | 6 | ≤ 3.24° vs. encoders | ≤ 4.27 Nm vs. torque sensors |
| KUKA iiwa Robot (IMU + markers) | 6 | ≤ 2.84° vs. encoders | ≤ 3.71 Nm vs. torque sensors |

## Measurement Configurations

The algorithm supports different sensor combinations via the measurement model:

- **TCMA_I** — IMU accelerometer and gyroscope only
- **TCMA_I&M** — IMU + optical motion capture markers
- **TCMA_I,M&0T** — IMU + markers + virtual zero-torque measurements

Additional modalities (e.g., force plates, markerless video) can be added by defining new measurement models in terms of the state vector.

## Tuning Parameters

- **Process noise `Q_τ`**: Controls how freely torques can change between time steps. Set as a diagonal matrix; values around 1 (Nm)² worked well across both experimental setups.
- **Measurement noise `R`**: Constructed from sensor residual covariances. The algorithm showed low sensitivity to precise values of `R`.
- **IEKF iterations (`NoI`)**: Number of measurement update iterations per time step (default: 3).

## Citation

If you use this code in your research, please cite:

```bibtex
@article{osman2026tcma,
  title={A Tightly Coupled Multibody Dynamics and Multi-Sensor Fusion Algorithm for Simultaneous Kinematics and Kinetics Estimation},
  author={Osman, Hassan and de Kanter, Daan and Boelens, Jelle and Kok, Manon and Seth, Ajay},
  year={2026}
}
```

## License

This project is licensed under the [Apache License 2.0](LICENSE).

## Acknowledgments

- Pendulum experimental data provided by Eyal Bar-Kochba, Johns Hopkins Applied Physics Laboratory.
- This work was supported by the Convergence Human Mobility Centre, funded by the Convergence Alliance (TU Delft, Erasmus MC, Erasmus University Rotterdam).
- Built on [OpenSim](https://opensim.stanford.edu/) and [Simbody](https://simtk.org/projects/simbody).
