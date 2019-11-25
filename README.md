# Spectral Neighbor Analysis Potential (SNAP)

This directory contains a spectral neighbor analysis potential (SNAP[1]) Driver. It uses bispectrum components to characterize the local neighborhood of each atom in a very general way.

## SNAP model driver

The SNAP driver is written in C++ and it expects a SNAP coefficient file followed by a SNAP parameter file.

### Coefficient file format

The format of the SNAP coefficient file is as follows:

Blank lines and lines beginning with the `#` character are ignored.

```bash
nelem, ncoeff
```

where, `nelem` is the number of elements, and `ncoeff` is the number of coefficients.

It follows with one block for each of the `nelem` elements. Where, the first line of each block contains:

```bash
Element, R, W
```

where, `Element` is the element symbol (text string), `R` is the element radius (in distance units), and `W` is the element weight (dimensionless).

Lines `2, 3,..., ncoeff+1` of each block is

```bash
coeffs
```

where, `coeffs` are SNAP coefficients, one per line.

### Parameter file format

The format of the SNAP parameter file is as follows:

Blank lines and lines beginning with the `#` character are ignored.\
Non-blank and non-comment lines must contain one keyword/value pair.

```bash
keyword value
```

The mandatory keywords are `rcutfac` and `twojmax`.
Optional keywords are `rfac0`, `rmin0`, `switchflag`, `bzeroflag`, and `quadraticflag`.

## Hybrid style

SNAP can be used in a hybrid style and in combination with a Ziegler-Biersack-Littmark (ZBL[2,3]) screened nuclear repulsion potential and any other two-body interaction in the table style[3,4] representation.

### HYBRID parameter file format

The format of the HYBRID parameter file is as follows:

Blank lines and lines beginning with the `#` character are ignored.

__NOTE:__ __`No mixing rule will be used here.`__

Each pair `(i,j)` or `(j,i)` can be assigned to one style. If you specify the same pair for the second time, it wipes
out all the previous assignments of that pair and the second one will be calculated for the two interacting atoms
of those types.

```bash
N
```

`N` is the number of elements in the hybrid style

```bash
Element_1 Element_2 ... Element_N
```

`Element_1 Element_2 ... Element_N` are number of elements names (atom names).

If there is any `zbl` interaction,

```bash
zbl  inner  outer
```

In the ZBL style, the inner and outer cutoff are the same for all pairs of atom types. It follows by,

```bash
Element_i  Element_j  zbl  Z_i  Z_j
```

where the `Element_i` and `Element_j` are the `i` and `j` element names respectively.
The values of `Z_i` and `Z_j` are equal to the atomic numbers of the two atom types.

If there is any `table` interaction,

```bash
table  style  N
```

where `style` is the method of interpolation and one of the `lookup` or `linear` or `spline` or `bitmap`.\
`N` means to use `N` values in either of the `lookup`, `linear`, `spline` table styles, or\
`N` means to use `2^N` values in the `bitmap` table style.

```bash
Element_i  Element_j  table  style_number  filename  keyword  [cutoff]
```

`table` style can be used multiple times. For example, if for the interactions between `i` and `i` atoms we use a `linear` table style and for the interactions between `i` and `j` atoms, we use a `spline` table style, then you should list the table style two times as below:

```bash
table  linear  1000
table  spline  2000
```

where, style indexing starts from `1`. It means, that the fisrt style is numbered as `1` and the second one as `2` and so on so forth. Later in the pair interactions, the `table` style must be added after the `i`, and `j` atom names followed by the style number and then followed by the remaining coefficients as of `filename`, `keyword`, and/or `cutoff`.\
The `filename` specifies a file containing tabulated energy and force values.\
The `keyword` specifies a section of the file.\
The `cutoff` is an optional coefficient in distance unit.

__NOTE:__ __`keyword` should be unique for each pair interaction__

Example `1`, where we use 2 styles

```bash
table  linear  1000
table  spline  10000
i  i  table  1  ii.table ii_keyword 4.0
i  j  table  2  ij.table ij_keyword
```

Example `2`, we only use one style

```bash
table  spline  10000
i  i  table  1  ii.table ii_keyword 4.0
i  j  table  1  ij.table ij_keyword
```

Example `3`, we use one style and only have one tabulated file

```bash
table  spline  10000
i  i  table  1  tablefile.txt ii_keyword 4.0
i  j  table  1  tablefile.txt ij_keyword 4.8
```

## References

1. [Thompson, Swiler, Trott, Foiles, and Tucker, "Spectral neighbor analysis
method for automated generation of quantum-accurate interatomic potentials,"
J Comp Phys, 285, 316 (2015)](https://www.sciencedirect.com/science/article/pii/S0021999114008353).

2. J.F. Ziegler, J. P. Biersack, and U. Littmark, “The Stopping and Range of
Ions in Matter,” Volume 1, Pergamon, 1985.

3. [https://lammps.sandia.gov](https://lammps.sandia.gov)

4. [https://lammps.sandia.gov/doc/pair_table.html](https://lammps.sandia.gov/doc/pair_table.html)

## Contributing

Copyright (c) 2019, Regents of the University of Minnesota.\
All rights reserved.

Contributors:\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Yaser Afshar

## License

[LGPLv2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
