# Flat plate transitional flows - T3A case

For reference see classical ERCOFTAC test cases at [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020)

The case is solved in a rectangular domain spanning from $x = -0.15$ m to $x = 1.7$ m. The height of the domain is $h - 0.3$ m and the leading edge of the plate is at $x=0$ m.

Three different meshes were prepared: coarse, medium, and fine. 

| Mesh | Size | average $y^+_1$ |
|------|------|-----------------|
| coarse | 157 x 125 |          |
| medium | 314 x 250 |          |
| fine   | 618 x 500 |          |

An incompressible fluid with constant density and constant viscosity $\nu=1.5\times 10^{-5}$ m2/s is assumed.

Boundary conditions are set as follows:
* inlet ($x=-0.15$ m)
    - $U=5.4\,m/s$
    - $Tu=3.3\,\%$
    - $\nu_t/\nu = 12$
    
One can run the simulation with 

```
./Allrun medium
```
