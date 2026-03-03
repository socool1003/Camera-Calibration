This project aims to implement and compare two classic camera calibration methods, **Zhang's Method** and **Tsai's Method**, using MATLAB. The project first generates simulated camera capture data, then applies both calibration methods to estimate the camera's intrinsic and extrinsic parameters, followed by result analysis and visualization.

## Project Overview

Camera calibration is a fundamental task in computer vision, used to determine a camera's intrinsic parameters (such as focal length, principal point, distortion coefficients) and extrinsic parameters (the camera's pose in the world coordinate system). The core logic of this project is divided into three main parts:

1.  **Synthetic Data Generation**: Creation of a synthetic dataset with "ground truth" parameters, serving as input and validation basis for both calibration algorithms.
2.  **Zhang's Calibration Method Implementation**: Following Zhang's classic four-step approach, from homography matrix calculation to non-linear optimization, to progressively estimate camera parameters.
3.  **Tsai's Calibration Method Implementation**: Utilizing Tsai's linear two-step method to efficiently estimate camera parameters from a single image.

By implementing these two classic methods, the project aims to deepen the understanding of their mathematical principles, computational procedures, and their respective characteristics and application scenarios.

---

## Part 1: Synthetic Data Generation (`gen_synthetic_data.m`)

This script is responsible for creating a controllable synthetic dataset for camera calibration algorithms. It simulates a camera capturing a standard chessboard pattern from various poses and calculates the projected pixel coordinates of each chessboard corner in the image plane.

**Core Logic:**
*   Defines the true camera intrinsic parameters `K_true` and a 3D chessboard model.
*   Randomly generates multiple camera poses (rotation `R` and translation `t`), simulating multi-view capture.
*   Calculates the precise mapping from world points to pixel points based on geometric projection principles.
*   Saves the generated `K_true`, `world_points`, `image_points`, and `extrinsics` (true extrinsic parameters) to `synthetic_calib_data.mat`, serving as input and validation benchmark for subsequent calibration steps.

---

## Part 2: Zhang's Calibration Method (Zhang's Method)

Zhang's calibration method uses multiple images of a planar chessboard to solve for camera parameters. This implementation is divided into four sub-steps, progressively recovering camera intrinsic and extrinsic parameters from image points:

1.  **Step 1: Calculate Homography Matrix H (`zhang_calibration_step1.m`)**
    *   **Purpose**: To calculate the $3 \times 3$ homography matrix `H` that connects the world coordinate system (chessboard plane) to the image pixel coordinate system for each image.
    *   **Method**: Utilizes Direct Linear Transform (DLT) to construct a system of linear equations and solves it via Singular Value Decomposition (SVD).
2.  **Step 2: Solve for Intrinsic Parameters K (`zhang_calibration_step2.m`)**
    *   **Purpose**: To calculate the camera's intrinsic parameter matrix `K` in a closed-form solution, leveraging the geometric constraints (orthogonality of the chessboard) embedded in the homography matrices.
    *   **Method**: Extracts constraints from `H`, constructs a linear system `V*b = 0`, solves for matrix `B`, and then inverts to find the elements of `K`.
3.  **Step 3: Solve for Extrinsic Parameters (R, t) and Visualize (`zhang_calibration_step3.m`)**
    *   **Purpose**: Given `K` and `H`, to decompose `H` to extract the camera's extrinsic parameters (rotation matrix `R` and translation vector `t`) corresponding to each image.
    *   **Method**: Uses `K_inv * H` to decompose `r1, r2, t`, computes `r3` via cross product, and applies SVD for orthogonalization correction. Visualizes the comparison between true and estimated camera poses in 3D.
4.  **Step 4: Non-linear Optimization (Bundle Adjustment) (`zhang_calibration_step4.m`)**
    *   **Purpose**: To globally refine all camera parameters (intrinsic, extrinsic, and distortion coefficients) by minimizing the reprojection error across all images.
    *   **Method**: Employs the Levenberg-Marquardt algorithm (`lsqnonlin`) for iterative optimization to find the optimal parameters and visualizes the reprojection error before and after optimization.

---

## Part 3: Tsai's Calibration Method (Tsai's Method) (`tsai_calibration.m`)

Tsai's calibration method is an efficient, linear two-step approach that can theoretically estimate camera intrinsic and extrinsic parameters from a single image.

**Core Logic:**
*   **Image Coordinate Preprocessing**: Converts raw pixel coordinates to image physical coordinates centered at the principal point.
*   **Step 1: Solve for R and tx, ty based on Radial Alignment Constraint (RAC)**: Utilizes the geometric property of the radial alignment constraint to establish a system of linear equations, and solves for partial elements of the rotation matrix `R` and the `x, y` components of the translation vector `t` using least squares.
*   **Step 2: Solve for tz and f**: Based on the results from Step 1, constructs another system of linear equations to linearly solve for the `z` component of the translation vector and the camera's focal length `f`.
*   **Result Comparison and Visualization**: Compares the estimated focal length and translation vector with their true values, and visualizes the camera pose recovery effect in 3D.

---

## Technical Stack & Dependencies

*   **Programming Language**: MATLAB
*   **Core Libraries**: MATLAB built-in functions (`svd`, `norm`, `cross`, `inv`, `lsqnonlin`, etc.).
*   **File Dependencies**: Scripts pass data via `.mat` files; ensure they are run in sequential order.

## How to Run

1.  **Environment Setup**: Ensure MATLAB is installed on your system, and the **Optimization Toolbox** is also installed (required for the non-linear optimization step in Zhang's Method).
2.  **File Placement**: All `.m` script files should be placed in the same working directory in MATLAB.
3.  **Execution Order**: Please run the scripts in the following exact order from the MATLAB command line:
    1.  `gen_synthetic_data`
    2.  `zhang_calibration_step1`
    3.  `zhang_calibration_step2`
    4.  `zhang_calibration_step3`
    5.  `zhang_calibration_step4`
    6.  `tsai_calibration`
4.  **View Output**: After each script runs, the MATLAB command window will display calculation results, error analysis, and pop up corresponding visualization plots.
