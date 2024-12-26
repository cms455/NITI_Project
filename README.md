# Supplementary Information
*By Calvin Smith*  
*December 2024*

## S1: Analytical model of corneal elasticity

Our goal is to characterize the Lamé parameters, G and μ, of the compliance matrix of the cornea. The values of G and μ will tell us about the transverse anisotropy of the cornea and help characterize ocular disease. Simulated dispersion curves from Pitre *et al.* are used to fit the data collected from the OCE systems using a simple least squares linear regression.

### Data Measurement

The OCE method used in this project utilizes only four frequency measurements of the Cornea to decrease measurement time and improve recovered signal. As a result our code takes as input data four frequencies in the frequency/wavenumber plane. For each frequency (600, 900, 1200, 1500 Hz) we get a range of wavenumbers. The points we use to fit the data are the peak k-values along the A0 mode. Our final input is 4 points representing the wavenumber peaks along four select frequencies.

### Physical Model

Following Pitre *et al.*, we model the cornea as a nearly-incompressible transverse isotropic (NITI) material. This model was validated in porcine retinas ex vivo by Pitre *et al.* Hooke's law for a NITI material can be expressed in the Voigt notation via the stiffness matrix C:

```math
\begin{bmatrix}
    \sigma_{xx}
    \\ \sigma_{yy}
    \\ \sigma_{zz}
    \\ \tau_{yz}
    \\ \tau_{xz}  
    \\ \tau_{xy}
\end{bmatrix} =
\begin{bmatrix}
    \lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\
    \lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\
    \lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\
    0 & 0 & 0 & G & 0 & 0 \\
    0 & 0 & 0 & 0 & G & 0 \\
    0 & 0 & 0 & 0 & 0 & \mu
\end{bmatrix}
\begin{bmatrix}
    \varepsilon_{xx}
    \\ \varepsilon_{yy}
    \\ \varepsilon_{zz}
    \\ \gamma_{yz}
    \\ \gamma_{xz}  
    \\ \gamma_{xy}
\end{bmatrix}
```

where σ, τ, γ, ε are the stress and strain tensors, λ and μ are Lamé's parameters and the only anisotropy resides in the shear term C₄₄ = G -- the material is assumed isotropic in its longitudinal terms. For the cornea, following Pitre *et al.*, we use λ = ρcₚ² - 2μ, where ρ = 1000 kg/m³ is the cornea density, and cₚ = 1540 m/s is the corneal speed of sound.

We utilize the following parameters as inputs into the Pitre Dispersion Curve Simulations:

- ρ = 1000: Density of the material (kg/m³)
- ρₗ = 1000: Density of the Liquid (Water) (kg/m³)
- cₗ = 1480: Speed of Sound in the Liquid (m/s)
- cₚ = 1540: Speed of Sound in the Tissue (m/s)
- h = 550 μm: Thickness (μm)

### Boundary Conditions

We assume a solid planar cornea of thickness h bounded by semi-infinite media with vacuum on top and water, modeled as an acoustic fluid (ρₗ = 1000 kg/m³, cₗ = 1480 m/s).

![Geometry and Boundary Conditions of the Cornea](./images/Cornea_diagram.png)

Dimensionless parameters are defined as:
```math
\begin{aligned}
x' &= x/h \\
u' &= u/h \\
t' &= \frac{t}{h} \sqrt{\frac{\mu}{\rho}} \\
f' &= t \cdot h \sqrt{\frac{\mu}{\rho}} \\
k' &= k \cdot h
\end{aligned}
```

The general solution for the dimensionless horizontal and vertical displacements (u,v) in the solid is a linear combinations of harmonic solutions described as:

```math
\begin{aligned}
u(x,z,t) &= \sum_{j=1}^4 C_jA_je^{il_{j}z}e^{i(kx-wt)} \\
v(x,z,t) &= \sum_{j=1}^4 C_je^{il_{j}z}e^{i(kx-wt)}
\end{aligned}
```

where:

```math
A = \pm \frac{\sqrt{2} \frac{\gamma^2 k}{\alpha^2} \sqrt{\phi \pm\sqrt{\phi^2 - 4q_{\alpha}^2q_{\beta}^2}}}{\phi + 2\frac{\beta^2}{\alpha^2}q_{\beta}^2\pm \sqrt{\phi^2 - 4q_{\alpha}^2q_{\beta}^2}}
```

```math
l = \pm \sqrt{\frac{1}{2}[\phi \pm \sqrt{\phi^2 - 4q_{\alpha}^2q_{\beta}^2}]}
```

```math
\begin{aligned}
q_{\alpha}^2 &= k^2 - \frac{w^2}{\alpha^2} \\
q_{\beta}^2 &= k^2 - \frac{w^2}{\beta^2} \\
\phi &= \frac{\gamma^4k^2}{\alpha^2\beta^2} - \frac{\beta^2}{\alpha^2}q_{\alpha}^2 - \frac{\beta^2}{\alpha^2}q_{\beta}^2
\end{aligned}
```

```math
\alpha^2 = \frac{G}{\mu}, \quad \beta^2 = \frac{\lambda + 2\mu}{\mu}, \quad \gamma^2 = \frac{\lambda + G}{\mu}
```

The fluid-domain wave equation is expressed in the velocity potential form:
```math
\begin{aligned}
\mathbf{\dot u}^f &= \nabla \Phi \\
p^f &= \frac{\rho^f}{\rho}\Phi_t \\
\Delta \Phi &- \frac{1}{\delta^2}\Phi_{tt} = 0 \\
\delta^2 &= \frac{\rho(c^f)^2}{\mu}
\end{aligned}
```

The general solution in the fluid domain is:
```math
\begin{aligned}
\Phi &= C_5e^{\xi z}e^{i(kx-wt)} \\
\xi &= \sqrt{k^2 - w^2/\delta^2} \\
\text{Re}(\xi) &> 0
\end{aligned}
```

where Φ is the velocity potential and δ² = ρ(cₗ)²/μ

The constants Cⱼ are found by substituting the general solutions into the boundary conditions:
- σₓₓ = 0
- σᵧᵧ = 0 at z = 1
- σₓᵧ = 0 at z = 0
- σᵧᵧ = σᵧᵧᶠ at z = 0
- v̇ = v̇ᶠ at z = 0

This can be rewritten as **M**c = 0, where c = [C₁, C₂, C₃, C₄, C₅]. This system of equations can be numerically solved by minimizing the absolute value of det **M**. This gives rise to symmetric and antisymmetric waveguide modes (similar to Lamb modes) that are characterized by their dispersion curves k(w).

### Least Squares Linear Regression

Given our data, which is only 4 points in the frequency-wavenumber space, we experimented and found that a simple linear regression model using least squares had the best results. Based on a starting G₀ and μ₀, we generate a dispersion curve using the model developed by the Pitre Group. Our loss function is defined as:

```math
\text{Loss} = \sum_{i=1}^4 ( \text{Data}(f_i)- D(f_i, G ,\mu))^2
```

Where D(f, G, μ) is the calculated dispersion curve, Data(f) is the data, and f⃗ = [f₁, f₂, f₃, f₄] are the selected frequencies. Then we utilize the Matlab tool fminsearch to minimize the calculate loss as function of G and μ.

<div style="display: flex; justify-content: space-between;">
  <div style="flex: 1;">
    <img src="./images/linear_fit_v1.jpg" alt="Linear Regression fit for simulated data with no noise" width="90%"/>
    <p><strong>Fig 1.</strong> Linear Regression fit for simulated data with no noise.<br/>
    Initial Guess: G₀ = 30E3, μ₀ = 30E3<br/>
    Fit D-Curve: G = 30E3, μ = 30E3</p>
  </div>
  <div style="flex: 1;">
    <img src="./images/linear_fit_w_noise_v1.jpg" alt="Linear Regression Fit for simulated data with added noise" width="90%"/>
    <p><strong>Fig 2.</strong> Linear Regression Fit for simulated data with added noise of SNR = 20<br/>
    Initial Guess: G₀ = 30E3, μ₀ = 30E3<br/>
    Fit D-Curve: G = 29.8E3, μ = 162.6E3</p>
  </div>
</div>

## Results

Looking at the dispersion curves for a range of G and μ values between [2E4, 2E6], which is a reasonable range for the Lamé Parameters in the Cornea, we see that the lines converge quickly as G and μ increase (Fig 3., Fig 4.).

To visualize the robustness of the linear fit to noise, we created a contour plot of the gaussian distribution of G and Mu values for SNR of 20. There is significant clustering around the actual G and Mu values, however there are significant numbers of outliers and a lack of precision. The reason for this is most likely due to random noise on the 4 data points averaging out to lowering the fit, causing the G and μ parameters to become unstable, with values rapidly increasing to unreasonable values.

![Dispersion curves of various G values](./images/GPlot_v4.jpg)
**Fig 3.** Dispersion curves of various G values for μ₀ = 4e6. The lines represent the dispersion curves of the G values between [20kpa, 2000kpa].

![Dispersion curves of various μ values](./images/Muplot_v2.jpg)
**Fig 4.** Dispersion curves of various μ values for G₀ = 25.6 kpa. The lines represent the dispersion curves of the μ values between [20kpa, 2000kpa].

![Fitted G and μ values](./images/banana_plot_v1.jpg)
**Fig 5.** Fitted G and μ values for G = 30E3, μ = 30E3 with SNR = 20 for 100 iterations.

