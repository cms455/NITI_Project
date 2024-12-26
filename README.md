# Supplementary Code for NITI OCE Linear Regression Fitting
*By Calvin Smith*  
*December 2024*

## S1: Analytical model of corneal elasticity

Our goal is to characterize the Lamé parameters, G and μ, of the compliance matrix of the cornea. The values of G and μ will tell us about the transverse anisotropy of the cornea and help characterize ocular disease. Simulated dispersion curves from Pitre *et al.* are used to fit the data collected from the OCE systems using a simple least squares linear regression.

### Data Measurement

The OCE method used in this project utilizes only four frequency measurements of the Cornea to decrease measurement time and improve recovered signal. As a result our code takes as input data four frequencies in the frequency/wavenumber plane. For each frequency (600, 900, 1200, 1500 Hz) we get a range of wavenumbers. The points we use to fit the data are the peak k-values along the A0 mode. Our final input is 4 points representing the wavenumber peaks along four select frequencies.

### Physical Model

We utilize the following parameters as inputs into the Pitre Dispersion Curve Simulations:

- ρ = 1000: Density of the material (kg/m³)
- ρₗ = 1000: Density of the Liquid (Water) (kg/m³)
- cₗ = 1480: Speed of Sound in the Liquid (m/s)
- cₚ = 1540: Speed of Sound in the Tissue (m/s)
- h = 550 μm: Thickness (μm)

The physics of the cornea is defined using the Nearly-Incompressible-Transverse-Isotropy (NITI) model. The NITI compliance matrix is calculated using two moduli (μ, G) where λ = ρcₚ² - 2μ and is defined as:

```math
\mathbf{C} =
\begin{bmatrix}
    \lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\
    \lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\
    \lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\
    0 & 0 & 0 & G & 0 & 0 \\
    0 & 0 & 0 & 0 & G & 0 \\
    0 & 0 & 0 & 0 & 0 & \mu
\end{bmatrix}
```

### Boundary Conditions

The final boundary condition in the fluid domain defines a general solution of:

- Φ = C₅e^(ξz)e^(i(kx-wt))
- ξ = √(k² - w²/δ²)
- Φ = velocity potential, δ² = ρ(cₗ)²/μ

Applying the BCs to the EOMs results in a linear set of equations described by a characteristic matrix **M** such that **M**c = 0, where c = [C₁, C₂, C₃, C₄, C₅]. Since **M**c = 0, by the invertibility of the determinant then det**M** = 0. The Pitre Group utilizes the explained method to calculate the dispersion curve for a given G and μ by finding the wave-number k for each frequency that minimizes |det**M**|.

### Least Squares Linear Regression

Given our data, which is only 4 points in the frequency-wavenumber space, we experimented and found that a simple linear regression model using least squares had the best results. Based on a starting G₀ and μ₀, we generate a dispersion curve using the model developed by the Pitre Group. Our loss function is defined as:

```math
\text{Loss} = \sum_{i=1}^4 ( \text{Data}(f_i)- D(f_i, G ,\mu))^2
```

Where D(f, G, μ) is the calculated dispersion curve, Data(f) is the data, and f = [f₁, f₂, f₃, f₄] are the selected frequencies. Then we utilize the Matlab tool fminsearch to minimize the calculate loss as function of G and μ.

## Results

Looking at the dispersion curves for a range of G and μ values between [2E4, 2E6], which is a reasonable range for the Lamé Parameters in the Cornea, we see that the lines converge quickly as G and μ increase (Fig 3., Fig 4.).

To visualize the robustness of the linear fit to noise, we created a contour plot of the gaussian distribution of G and Mu values for SNR of 20. There is significant clustering around the actual G and Mu values, however there are significant numbers of outliers and a lack of precision. The reason for this is most likely due to random noise on the 4 data points averaging out to lowering the fit, causing the G and μ parameters to become unstable, with values rapidly increasing to unreasonable values.

### Key Findings:

1. Linear regression fit for simulated data with no noise:
   - Initial Guess: G₀ = 30E3, μ₀ = 30E3
   - Fit D-Curve: G = 30E3, μ = 30E3

2. Linear regression fit for simulated data with added noise (SNR = 20):
   - Initial Guess: G₀ = 30E3, μ₀ = 30E3
   - Fit D-Curve: G = 29.8E3, μ = 162.6E3


