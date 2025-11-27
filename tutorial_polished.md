# Tutorial: Reduced-Order Flexible Solid Modelling for Rocking Problems in Simulink

## Overview
This tutorial provides the MATLAB codes, the Simulink models, and the necessary introductions required to run the examples presented in the companion paper. Two examples from the companion paper are covered:

1. A realistic column-like structure  
2. A scaled frame structure used in a physical experiment  

The first example demonstrates how to use MATLAB and Simulink to conduct reduced-order analysis for rocking problems.  
The second example shows how to build a structural system with many inter-connected components in Simulink.

---

## Pre-requisites

Readers must have a good understanding of:

1. The **Craig-Bampton reduction method**.  
2. The **Floating Frame of Reference formulation**.

### Suggested Readings
- The companion paper, which provides an introduction to both theories.  
- The paper **“[Coupling of Substructures for Dynamic Analyses](https://hal.science/hal-01537654v1/document)”**, which orginally proposed the Craig-Bampton method.  
- The paper **“[Primer on the Craig-Bampton](https://www.vibrationdata.com/tutorials2/Primer_on_the_Craig-Bampton_Method.pdf)”**, which provides a well-written introduction to the Craig-Bampton method, including checks to verify successful reduction and other practical tips.  
- **Chapter 5** in the book **“[Dynamics of Multibody Systems](https://www.cambridge.org/highereducation/books/dynamics-of-multibody-systems/E287DA737B6138E040AA96FC12F7D7DF#contents)”**, which proposes and elaborates on the floating frame of reference formulation.
- This well-written **[Simscape Example](https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html)**, which can make you familiar with the basic workflow and common functions for using the reduced order flexible solid in Simscape.

📌*I highly recommend reviewing the theory-related materials first to build a clear conceptual understanding before moving on to the Simscape Example. Otherwise, it may be difficult to distinguish between syntax-related logic and theory-related logic, which can lead to confusion.

### Software Requirements:

1. MATLAB 2023a or later version.
2. Simulink 2023a or later version.
3. MATLAB **[Partial Differential Equation Toolbox](https://uk.mathworks.com/products/pde.html)**.
4. Simulink extension **[Simscape Multibody](https://uk.mathworks.com/products/simscape-multibody.html)**.
5. OxContact library in Simscape. 
---



## Why This Separate Documentation?
Although MATLAB provides the **[Simscape Example](https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html)** for using the reduced order flexible solid, this separate documentation is necessary for the following reasons.

### 1. Additional functionalities not included in the Simscape Example
The provided codes add two important capabilities to MATLAB’s original implementation:

1. **Allows retention of any desired fixed-interface modes**, whereas MATLAB’s implementation only keeps fixed-interface modes within one specific frequency range.  
2. **Allows reconstruction of full-order deformation, strain, stress, and global displacements** (add credits to houmain, inform Manolis on this, GPT do not change this line) from Simulink simulation results, whereas MATLAB’s implementation only returns global displacements of the interface nodes in Simulink.

_The first example in this tutorial covers these functionalities._

---

### 2. The Simscape Example does not address rocking problems
The Simscape Example is not a rocking problem and therefore does **not** provide guidance on how to use the **OxContact elements** (an in-house toolbox for contact modelling in rocking problems)  

_Both examples address this._

---

### 3. The Simscape Example assumes rigid interfaces
The Simscape Example applies only to structural components where the interface is stiff and can safely be assumed to be rigid. However, rocking components **without** interface reinforcement cannot be assumed to have **rigid interfaces**, requiring additional care.

_The first example exclusively covers this scenario._

---

### 4. The Simscape Example involves only a few physical components
Civil engineering structures such as frame structures usually include many inter-connected components. Constructing such a system in Simulink — which uses the floating frame of reference formulation — requires more effort than traditional civil engineering software.

_The second example in this tutorial provides a reference case._

---

### 5. Craig-Bampton and FFR are relatively unfamiliar in civil engineering
Despite their advantages for addressing rocking problems, the:

- **Floating frame of reference formulation**, and  
- **Craig-Bampton reduction method**

remain relatively unfamiliar to many civil engineerin engineers and researchers.

_The examples in this tutorial aim to improve familiarity with these modelling techniques._

---

## General Workflow

The overall workflow is as follows:

1. **Perform Craig-Bampton reduction in MATLAB**  
2. **Extract reduced-order stiffness, mass, and damping matrices**  
3. **Transfer the reduced-order data to the “Reduced-Order Flexible Solid” block in Simulink**, which uses the floating frame of reference formulation  
4. **Run dynamic simulations in Simulink**  
5. **Reconstruct full-order solutions** from the reduced simulation results  

---

## Example 1: Rocking Column with Corner Feet
[Link text](./Example1.md)
