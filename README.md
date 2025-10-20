# Batch Gradient Descent (BGD) Unitary Optimisation — Octave Implementation

This repository contains an **Octave implementation** of a **Batch Gradient Descent (BGD)** algorithm designed to approximate a *reference unitary matrix* `Uref` through numerical optimisation.  
The method uses finite-difference gradients to minimise the Frobenius norm squared difference between `Uref` and a trial unitary matrix `Utest = UC(P)`.

The code produces detailed convergence plots, percentile analyses, and saves both the numerical results and the best unitary matrix obtained during optimisation.

---

## 🧠 Overview of the Algorithm

The optimisation follows these main steps:

1. **Reference Unitary Generation**  
   A random complex matrix `A` is generated and decomposed using a QR factorisation:
   \[
   U_\text{ref} = Q \cdot \frac{\text{diag}(R)}{|\text{diag}(R)|}
   \]
   ensuring that `Uref` is unitary.

2. **Optimisation Loop (Gradient Descent)**  
   For each simulation:
   - Initialise a random real matrix `P`.
   - Map `P` to a unitary matrix `Utest = UC(P)` via QR decomposition.
   - Evaluate the cost function:
     \[
     C = ||U_\text{ref} - U_\text{test}||_F^2
     \]
   - Compute the gradient of the cost function using **finite differences**:
     \[
     \frac{\partial C}{\partial P_{ij}} = 
     \frac{C(P_{ij} + \varepsilon) - C(P_{ij} - \varepsilon)}{2 \varepsilon}
     \]
   - Update the matrix `P` via gradient descent:
     \[
     P = P - \alpha \, \nabla_P C
     \]

3. **Convergence Tracking**
   - The algorithm records the cost at each iteration.
   - Convergence is reached when the cost drops below `0.01`.
   - Cost histories, percentiles, and iteration counts are plotted for visual inspection.

4. **Outputs**
   - Average and percentile cost plots.
   - Iteration counts until convergence.
   - The best unitary matrix found (`best_U_overall`) and the corresponding cost.

---

## ⚙️ Parameters and Settings

| Parameter | Symbol | Default | Description |
|------------|--------|----------|--------------|
| Learning rate | α | `0.35` | Step size for gradient descent updates. Affects convergence speed and stability. |
| Epsilon | ε | `1e-6` | Perturbation step used in finite-difference gradient calculation. |
| Iterations | *num_iter* | `50` | Maximum number of gradient descent iterations per simulation. |
| Random seed | — | `1111` | Fixed seed for reproducibility (`rng(1111)`). |
| Simulations | *num_simulations* | user-defined | Number of independent optimisation runs. |

---

### 🧮 Notes on Learning Rate (α)

- The default `α = 0.35` works well for most cases (`d = 2–5`).  
- **Reducing α** (e.g. to `0.25`) can improve stability for higher-dimensional matrices.  
- **Increasing α** may accelerate convergence but risks overshooting minima.  
- The learning rate is easy to change:
  ```matlab. python.. etc.
  learning_rate = 0.35;  % modify as needed

 ## 🧩 External Code Attribution

This repository includes the `UC.m` function authored by **Christoph Spengler (University of Vienna, 2011)**,  
which implements the *Composite Parameterisation of the Unitary Group U(d)*.

The function is provided under its original academic copyright and  
is used here solely for research and educational purposes.

If you use or reference this repository in your work, please cite the following papers as requested by the author:

> **Spengler, C., Huber, M., & Hiesmayr, B. C.** (2010).  
> *A composite parameterization of unitary groups, density matrices and subspaces.*  
> *Journal of Physics A: Mathematical and Theoretical,* 43 (38), 385306.  
> arXiv: 1004.5252  
>
> **Spengler, C., Huber, M., & Hiesmayr, B. C.** (2011).  
> *Composite parameterization and Haar measure for all unitary and special unitary groups.*  
> arXiv: 1103.3408

**Note:** The file `UC.m` remains © Christoph Spengler (2011) and is not distributed under the MIT License.  
All other code and documentation in this repository are © Liam (2025) and released under the **MIT License**.


