## Example 1: Rocking Column with Corner Feet

The first example models a **solid marble column** with four rigid corner feet standing on a high-grip, stiff support medium. The height of the corner feet is negligible.
This example also covers the detailed workflow, and the contents that are *not* covered by the **[Simscape Example](https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html)** start with ⭐ in this document.

Figure below shows the geometry of the column  
*(⚠️insert the figure)*

---

## Step 1: Craig-Bampton Reduction in MATLAB

### I — Define structure parameters and common properties  
(I will insert matlab code here)

---

### II — Generate or import the geometry  
Details can be found in the **[Simscape Example](https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html)**.

(I will insert matlab code here, do not change this line GPT)

---

### III — Mesh the component  
This step determines the position of the interface nodes.

Workflow:
- Mesh once to identify the nearby vertex locations  
- Add the necessary vertices  
- Re-mesh  
- Apply vertex-type multipoint constraints  

The multipoint constraint is required because tetrahedral element nodes have only **three DOFs**, while an interface in Simulink must have **six DOFs**. The multipoint constraint forms a small rigid plane that provides 6 DOFs at the interface.

(I will insert matlab code here, do not change this line GPT)

---

### ⭐IV — Locate nodes within the area of the feet  
(I will insert matlab code here, do not change this line GPT)

---

### ⭐V — Add the node positions as geometry vertices and re-mesh 
(I will insert matlab code here, do not change this line GPT)

---

### ⭐VI — Apply vertex-type multipoint constraints for the interface DOFs  
These DOFs are retained during Craig-Bampton reduction.

(I will insert matlab code here, do not change this line GPT)

---

### VII — Check the fixed-interface modes  
Take note of the frequencies corresponding to the modes you want to keep. To know what modes to keep, you need to run a convergence analysis with different modes until the results does not change with more modes. 

⭐Optionally, run the modal analysis within a **live script** and use the **[visualize PDE results task](https://uk.mathworks.com/help/pde/ug/visualizepderesults.html)** to visualize the mode shapes.

(I will insert matlab code here, do not change this line GPT)

*(⚠️insert the modeshape figure here)*

---

### VIII — Apply Craig-Bampton reduction  
⭐List the frequency ranges corresponding to the modes you want to retain.

(I will insert matlab code here, do not change this line GPT)

---

### IX — Manually calculate the reduced-order damping matrix  
This allows assignment of different modal damping ratios to different modes.

Optional:  
Simulink’s reduced-order flexible solid can automatically compute:

- uniform modal damping  
- Rayleigh damping  

(I will insert matlab code here, do not change this line GPT)

---
