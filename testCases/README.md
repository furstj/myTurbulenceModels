# Test Cases for myTurbulenceModels

**Note:** Tests are intended for validation/verification of turbulence models. The relaxation factors, tolerances, and schemes are not set to optimal values! 

## Flat plate transitional flows 

Classical ERCOFTAC test cases, see [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020)

- [T3A](Coupland/T3A/): $U=5.4$ m/s,  $Tu=3$ %, no pressure gradient
- [T3B](Coupland/T3B/): $U=9.4$ m/s,  $Tu=6$ %, no pressure gradient
- [T3A-](Coupland/T3AM/): $U=19.8$ m/s,  $Tu=0.9$ %, no pressure gradient
- [T3C1](Coupland/T3C1/): $U=5.9$ m/s,  $Tu=6.6$ %, with pressure gradient
- [T3C2](Coupland/T3C2/): $U=5.0$ m/s,  $Tu=3$ %, with pressure gradient
- [T3C3](Coupland/T3C3/): $U=3.7$ m/s,  $Tu=3$ %, with pressure gradient
- [T3C4](Coupland/T3C4/): $U=1.2$ m/s,  $Tu=3$ %, with pressure gradient
- [T3C5](Coupland/T3C5/): $U=8.4$ m/s,  $Tu=3$ %, with pressure gradient

* Flat plate transitional flows with separation

## Flat plate transitional flows with separation
see **Volino, R. J., & Hultgren, L. S.** (2001). *Measurements in Separated and Transitional Boundary Layers Under Low-Pressure Turbine Airfoil Conditions.* Journal of Turbomachinery, 123(2), 189–197. [https://doi.org/10.1115/1.1350408](https://doi.org/10.1115/1.1350408)

**Note:** Experimental data are obtained by digitization of graphs from literature. It is not accurate 

- [Re300k_lowTI](HulgrenVolino/Re300k_lowTI): $Re=300000$, $Tu=0.2$ %
- [Re300k_highTI](HulgrenVolino/Re300k_highTI): $Re=300000$, $Tu=2.5$ %
- [Re50k_lowTI](HulgrenVolino/Re50k_lowTI): $Re=50000$, $Tu=0.2$ %
- [Re50k_lowTI](HulgrenVolino/Re50k_highTI): $Re=50000$, $Tu=2.5$ %
