# Spectral Neighbor Analysis Potential (SNAP)

This directory contains a spectral neighbor analysis potential (SNAP[1, 5]) Driver. It uses bispectrum components to characterize the local neighborhood of each atom in a very general way.

## SNAP model driver

The SNAP driver is written in C++ and it expects a SNAP coefficient file followed by a SNAP parameter file.

1) The SNAP coefficient file, usually ends in the `.snapcoeff` extension.
2) The SNAP parameter file, usually ends in the `.snapparam` extension.

SNAP can be used in a hybrid style and in combination with a Ziegler-Biersack-Littmark (ZBL[2,3]) screened nuclear repulsion potential and any other two-body interaction in the table style[3,4] representation.

3) Hybrid parameter file, usually ends in the `.hybridparam` extension.
   It includes a Ziegler-Biersack-Littmark (ZBL) pair styles parameter, as well as pair interactions through table style(s).
4) The interpolation table file(s) from potential energy and force values listed in (a) file(s) as a function of distance.

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

For example, below, see the SNAP coefficient file for a [SNAP potential for Tantalum (Ta)](https://openkim.org/id/SNAP_ThompsonSwilerTrott_2015_Ta__MO_359768485367_000)
```bash
# DATE: 2014-09-05 CONTRIBUTOR: Aidan Thompson athomps@sandia.gov CITATION: Thompson, Swiler, Trott, Foiles and Tucker, arxiv.org, 1409.3880 (2014)

# LAMMPS SNAP coefficients for Ta_Cand06A

1 31
Ta 0.5 1
-2.92477
-0.01137
-0.00775
-0.04907
-0.15047
0.09157
0.05590
0.05785
-0.11615
-0.17122
-0.10583
0.03941
-0.11284
0.03939
-0.07331
-0.06582
-0.09341
-0.10587
-0.15497
0.04820
0.00205
0.00060
-0.04898
-0.05084
-0.03371
-0.01441
-0.01501
-0.00599
-0.06373
0.03965
0.01072
```

### Parameter file format

The format of the SNAP parameter file is as follows:

Blank lines and lines beginning with the `#` character are ignored.\
Non-blank and non-comment lines must contain one keyword/value pair.

```bash
keyword value
```

The mandatory keywords are:
- `rcutfac`
- `twojmax`

Optional keywords are:
- `rfac0`
- `rmin0`
- `switchflag`
- `bzeroflag`
- `quadraticflag`
- `chemflag`
- `bnormflag`
- `wselfallflag`

For example, below, see the SNAP parameter file for a [SNAP potential for Tantalum (Ta)](https://openkim.org/id/SNAP_ThompsonSwilerTrott_2015_Ta__MO_359768485367_000)
```bash
# DATE: 2014-09-05 CONTRIBUTOR: Aidan Thompson athomps@sandia.gov CITATION: Thompson, Swiler, Trott, Foiles and Tucker, arxiv.org, 1409.3880 (2014)

# LAMMPS SNAP parameters for Ta_Cand06A

# required
rcutfac 4.67637
twojmax 6

# optional

rfac0 0.99363
rmin0 0
bzeroflag 0
quadraticflag 0
```

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


Example `4`, below, see the SNAP HYBRID parameter file for a [SNAP potential for Tantalum (Ta)](https://openkim.org/id/SNAP_ThompsonSwilerTrott_2015_Ta__MO_359768485367_000)
```bash
# DATE: 2014-09-05 CONTRIBUTOR: Aidan Thompson athomps@sandia.gov CITATION: Thompson, Swiler, Trott, Foiles and Tucker, arxiv.org, 1409.3880 (2014)

# Definition of SNAP potential Ta_Cand06A
# Assumes 1 LAMMPS atom type


# Number of elements for the hybrid style
1

# Number of elements names, (atom names)
Ta

# zbl  inner  outer
zbl  4.0  4.8

# Element_1  Element_2  zbl  Z_1  Z_2
Ta  Ta  zbl  73  73
```

## References

1. [Thompson, A.P., Swiler, L.P., Trott, C.R., Foiles, S.M., and Tucker, G.J., "Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials," J Comp Phys, 285, 316 (2015)](https://www.sciencedirect.com/science/article/pii/S0021999114008353).

2. Ziegler, J.F., Biersack, J.P., and Littmark, U., "The Stopping and Range of
Ions in Matter," Volume 1, Pergamon, (1985)

3. [https://lammps.sandia.gov](https://lammps.sandia.gov)

4. [https://lammps.sandia.gov/doc/pair_table.html](https://lammps.sandia.gov/doc/pair_table.html)

5. [Cusentino, M.A., Wood, M.A., Thompson, A.P., "Explicit Multielement Extension of the Spectral Neighbor Analysis Potential for Chemically Complex Systems," J. Phys. Chem. A, 124, 5456 (2020)](https://pubs.acs.org/doi/full/10.1021/acs.jpca.0c02450)

## Contributing

Copyright (c) 2019--2020, Regents of the University of Minnesota.\
All rights reserved.

Contributors:\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Yaser Afshar

## License

[LGPLv2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
