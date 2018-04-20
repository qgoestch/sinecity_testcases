## sinecity_testcases

This contains all the python scritpts that have been used for the calculations of the
five test cases presented in the CFA2018 proceeding paper (see ressources folder).

The five test cases are 2D geometries related to outdoor sound propagation scenes.
The **grid convergence study** used for this **code verification** gives the 
**convergence rates** and the **observed orders of accuracy** for each tested numerical method.

The updates for the Finite Difference Time Domain (FDTD), Transmission Line Matrix (TLM) 
and Finite Difference Helmholtz (FDH) are available in the num_methods folder.


## Installation

Python 2.7.* together with `numpy`, `scipy`, `os`, `site` and `matplotlib` packages are 
required for this code.
All folders, except ressources, need to be downloaded and placed on the same level.

## Examples

Results and code architecture is described in the pdf forms placed in the ressources folder.
All cases are described and the solutions are provided.

## Tests

From any of the case*_*.py, figures from the CFA conference proceeding can be recovered and 
recalculated.

## Documentation

The documentation is located in the ressources foler. A html version is available, where 
all modules are described from the python docstrings placed within each function.

## License

`sinecity_testcases` is distributed under the BSD 3-clause license. See LICENSE for more information.

## Contributing

Contributors are very welcome.
