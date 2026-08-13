# Flat-plate transitional flow - T3C3 case

For reference, see the classical ERCOFTAC test case in the [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020).

The case is solved in a domain spanning from $x = -0.15$ m to $x = 1.7$ m. The leading edge of the plate is located at $x = 0$ m ant he height of the domain is taken from Suluksna et al. (2009) as

$$
h(y) = 0.3 \min\left(
      1.231*x^6 - 6.705*x^5 + 14.061*x^4 - 14.113*x^3 + 7.109*x^2 - 1.900*x + 0.950
      , 1
\right)
$$

Three different meshes were prepared: coarse, medium, and fine. Values of $y^+_1$ are calculated witt the `kOmegaSSTLM` turbulence model. The mesh details are summarized in the following table:

| Mesh | Size | min $y^+_1$ | max $y^+_1$ | avg $y^+_1$ |
|------|------|-------------|-------------|-------------|
| coarse | 157 x 125 |        |        |        |
| medium | 314 x 250 | $0.02$ | $0.41$ | $0.08$ |
| fine   | 618 x 500 |        |        |        |

An incompressible fluid with constant density and constant viscosity, $\nu = 1.5\times 10^{-5}$ m$^2$/s, is assumed.

The inlet boundary conditions at $x = -0.15$ m are as follows:

- $U=3.7$ m/s
- $Tu=3.1$ %
- $\nu_t/\nu = 6$


For example, the simulation can be run on the medium mesh with

```
./Allrun medium
```

## References

1. [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020)
2. **Suluksna, K., Dechaumphai, P., & Juntasaro, E.** (2009). *Correlations for modeling transitional boundary layers under influences of freestream turbulence and pressure gradient*. International Journal of Heat and Fluid Flow, 30(1), 66–75. [https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004](https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004)

