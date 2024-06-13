# Snowflake

[![DOI](https://zenodo.org/badge/792208776.svg)](https://zenodo.org/doi/10.5281/zenodo.11070766)

Evolution of twist-three parton distribution functions.

## Compilation and usage

The default compilator is *f95*, with *openmp* library. It could be modified in the first lines of the makefile.

For the compilation of the library do **make** .

The result of compilation is stored in /mod and /obj directories. To link with your code, compile with 

```
(... your compilation line ...) files.obj -I/Snowflake/mod
```

A single-file code can be compiled with Snowflake by typing
```
make program TARGET=your/code/file
```

For example
```
make program TARGET=prog/EXAMPLE.f90
/.a.out
```

will compile EXAMPLE.f90.

## Further information

The details of algorithm and implementation can be found in the publication e-Print: 2405.01162 (https://inspirehep.net/literature/2782890). Please, do not forget to cite it, if you use Snowflake. Some extra infromation about Snoflake can be found in manual/manual.pdf. For further questions and suggestions, please, contact me.
