# Flat-plate transitional flow - T3A case

For reference, see the classical ERCOFTAC test case in the [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020).

The case is solved in a rectangular domain spanning from $x = -0.15$ m to $x = 1.7$ m. The height of the domain is $h = 0.3$ m, and the leading edge of the plate is located at $x = 0$ m.

Three different meshes were prepared: coarse, medium, and fine. Values of $y^+_1$ are calculated witt the `kOmegaSSTLM` turbulence model. The mesh details are summarized in the following table:

| Mesh | Size | min $y^+_1$ | max $y^+_1$ | avg $y^+_1$ |
|------|------|-------------|-------------|-------------|
| coarse | 157 x 125 | $0.17$ | $0.94$ | $0.24$ |
| medium | 314 x 250 | $0.09$ | $0.56$ | $0.13$ |
| fine   | 618 x 500 | $0.06$ | $0.33$ | $0.07$ |

An incompressible fluid with constant density and constant viscosity, $\nu = 1.5\times 10^{-5}$ m$^2$/s, is assumed.

The inlet boundary conditions at $x = -0.15$ m are as follows:

- $U=5.4$ m/s
- $Tu=3.3$ %
- $\nu_t/\nu = 12$

The other boundary conditions are set as follows:

- flat plate ($x>0$ m, $y=0$ m): no-slip wall
- lower symmetry ($x<0$ m, $y=0$ m): symmetry 
- top wall ($y=0.3$ m): slip
- outlet ($x=1.7$ m): pressure outlet with $p=0$ Pa

The simulation can be run on the medium mesh with

```
./Allrun medium
```
