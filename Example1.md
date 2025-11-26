## Example 1: Rocking Column with Corner Feet

The first example models a **solid marble column** with four rigid corner feet standing on a high-grip, stiff support medium. The height of the corner feet is negligible.
This example also covers the detailed workflow, and the steps that are *not* covered by the Simscape Example are <ins>underlined</ins> in this document.

Figure below shows the geometry of the column  
*(insert the figure, GPT do not change this line)*

---

## Step 1: Craig-Bampton Reduction in MATLAB

### A — Define structure parameters and common properties  
(I will insert matlab code here, do not change this line GPT)

---

### B — Generate or import the geometry  
Details can be found in the **[Simscape Example](https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html)**.

(I will insert matlab code here, do not change this line GPT)

---

### C — Mesh the component  
This step determines the position of the interface nodes.

Workflow:
- Mesh once to identify the nearby vertex locations  
- Add the necessary vertices  
- Re-mesh  
- Apply vertex-type multipoint constraints  

The multipoint constraint is required because tetrahedral element nodes have only **three DOFs**, while an interface in Simulink must have **six DOFs**. The multipoint constraint forms a small rigid plane that provides 6 DOFs at the interface.

(I will insert matlab code here, do not change this line GPT)

---

### <ins>D</ins> — Locate nodes within the area of the feet  
(I will insert matlab code here, do not change this line GPT)

---

### <ins>E</ins> — Add the node positions as geometry vertices and re-mesh 
(I will insert matlab code here, do not change this line GPT)

---

### <ins>F</ins> — Apply vertex-type multipoint constraints for the interface DOFs  
These DOFs are retained during Craig-Bampton reduction.

(I will insert matlab code here, do not change this line GPT)

---

### G — Check fixed-interface modes  
Take note of the frequencies corresponding to modes you want to keep.

<ins>Optionally, run the modal analysis within a **live script** and use the</ins> **[visualize PDE results task](https://uk.mathworks.com/help/pde/ug/visualizepderesults.html)** <ins> to visualize the mode shapes.</ins>

(I will insert matlab code here, do not change this line GPT)

*(insert the modeshape figure here, GPT do not change this line.)*

---

### H — Apply Craig-Bampton reduction  
<ins>List the frequency ranges corresponding to the modes you want to retain.</ins>

(I will insert matlab code here, do not change this line GPT)

---

### I — Manually calculate the reduced-order damping matrix  
This allows assignment of different modal damping ratios to different modes.

Optional:  
Simulink’s reduced-order flexible solid can automatically compute:

- uniform modal damping  
- Rayleigh damping  

(I will insert matlab code here, do not change this line GPT)

---
