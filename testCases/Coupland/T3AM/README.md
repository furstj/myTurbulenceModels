# Flat-plate transitional flow - T3A- case

For reference, see the classical ERCOFTAC test case in the [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020).

The case is solved in a rectangular domain spanning from $x = -0.15$ m to $x = 2.9$ m. The height of the domain is $h = 0.3$ m, and the leading edge of the plate is located at $x = 0$ m.

Three different meshes were prepared: coarse, medium, and fine. Values of $y^+_1$ are calculated witt the `kOmegaSSTLM` turbulence model. The mesh details are summarized in the following table:

| Mesh | Size | min $y^+_1$ | max $y^+_1$ | avg $y^+_1$ |
|------|------|-------------|-------------|-------------|
| coarse | 157 x 125 | $0.36$ | $2.22$ | $0.58$ |
| medium | 314 x 250 | $0.18$ | $1.31$ | $0.29$ |
| fine   | 618 x 500 | $0.09$ | $0.78$ | $0.15$ |

An incompressible fluid with constant density and constant viscosity, $\nu = 1.5\times 10^{-5}$ m$^2$/s, is assumed.

The inlet boundary conditions at $x = -0.15$ m are as follows:

- $U=19.8$ m/s
- $Tu=0.874$ %
- $\nu_t/\nu = 8.72$

The other boundary conditions are set as follows:

- flat plate ($x>0$ m, $y=0$ m): no-slip wall
- lower symmetry ($x<0$ m, $y=0$ m): symmetry 
- top wall ($y=0.3$ m): slip
- outlet ($x=2.9$ m): pressure outlet with $p=0$ Pa

The simulation can be run on the medium mesh with

```
./Allrun medium
```

## References

1. [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020)
2. **Suluksna, K., Dechaumphai, P., & Juntasaro, E.** (2009). *Correlations for modeling transitional boundary layers under influences of freestream turbulence and pressure gradient*. International Journal of Heat and Fluid Flow, 30(1), 66–75. [https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004](https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004)
