# Helmsdale Geothermal Feasibility Study

## 2D Variable-Coefficient Advection–Diffusion Thermal Model

This repository contains a MATLAB implementation of a two-dimensional variable-coefficient advection–diffusion thermal model developed for the EARTH5016 Numerical Geodynamics assessment.

The model simulates conductive heat transport with internal radiogenic heat production and simplified vertical advective heat transport within a fault zone.

---

## 1. Repository Structure

### main.m

Primary finite-difference solver.

- Runs a 150 kyr transient simulation  
- Computes a quasi-steady thermal structure  
- Produces:
  - 2D temperature distribution with selected economic isotherms
  - 1D geothermal gradient profiles (Granite, Basin, Fault)
  - Comparison with He1 borehole observations

### sensitivity_study.m

Parameter study script.

- Executes four geological scenarios
- Varies:
  - Granite radiogenic heat production (Qr)
  - Fault permeability proxy (KD)
- Extracts the depth of the 100 °C isotherm at the fault

### README.md

Project documentation and execution instructions.

---

## 2. How to Run the Code

Both scripts are fully self-contained.  
No external data files are required.

### 2.1 Requirements

- MATLAB (tested on R2023 or later)
- No additional toolboxes required

### 2.2 Run the Baseline Model

1. Open MATLAB.
2. Set the Current Folder to this repository directory.
3. Open main.m.
4. Click Run.

Expected output:

- 2D temperature field with selected isotherms
- 1D geothermal gradient comparison and borehole validation

Execution time: approximately 1–2 minutes.

### 2.3 Run the Sensitivity Study

1. Open sensitivity_study.m.
2. Click Run.

Expected output:

- Comparison of the 100 °C isotherm position across scenarios
- Printed 100 °C depth values in the Command Window

Execution time: several minutes.

---

## 3. Governing Equation

rho Cp (dT/dt) = div(k grad T) − rho_f Cp_f v · grad T + Qr

Material properties vary spatially across granite, sedimentary basin, and fault domains.

---

## 4. Numerical Implementation

### 4.1 Spatial Discretization

- Grid: 141 × 61 nodes  
- Domain size: 14 km × 6 km  
- Method: Finite Difference Method

The conductive term is implemented in flux form using face-averaged thermal conductivities.

### 4.2 Advection

Advection is discretized using a first-order upwind scheme in the vertical direction.

Fluid velocity is parameterized as a permeability-scaled vertical proxy within the fault zone.

### 4.3 Time Integration and Stability

- Explicit forward Euler time integration  
- Total simulation time: 150 kyr  
- Time step determined using a CFL-based stability constraint  

---

## 5. Boundary Conditions

- Surface: Fixed temperature (10 °C)  
- Base: Prescribed conductive gradient (35 °C/km)  
- Lateral boundaries: Zero horizontal heat flux  

---

## 6. Parameter Study Design

The sensitivity analysis evaluates:

1. Granite radiogenic heat production (Qr)  
2. Fault permeability proxy (KD)

For each scenario, the depth of the 100 °C isotherm at the fault is extracted.

---

## 7. Reproducibility

- All material properties are defined directly in the scripts.
- No external input/output files are required.
- Results are deterministic for a given MATLAB version and grid resolution.

---

## 8. Intended Use

Submitted as part of the EARTH5016 individual consultancy report.

The model represents a simplified thermal framework suitable for feasibility-level geothermal assessment.
