[![DOI](https://zenodo.org/badge/781862015.svg)](https://doi.org/10.5281/zenodo.14825532)
![MIT License](https://img.shields.io/badge/license-MIT-blue)
# Martini 3 model of Graphene
<p align="center">
<img src="GRA_AA_CG_model.png" width="600">
</p>

Generates a Martini 3 model of both finite and infinite graphene sheet for running the molecular dynamics simulation with the Gromacs simulation package. The script outputs both the structure file (.gro) and a topology file (.itp).

## Requirements
Python 3 with `numpy` and `MDAnalysis`. The runtime dependencies are listed in `requirements.txt`.

## Usage
For generating the structure and topology of a finite graphene sheet, run
 ```python
 python Non-Periodic/martini3-graphene-topology.py -x [Length of the graphene sheet along x in nm] \
                                         -y [Length of the graphene sheet along y in nm] \
                                         -z [Length of the box along z in nm]
                                         -o [Name of the output for both structure and topology] \
                                         [--output-dir DIR] \
                                         [--force] \
                                         [--quiet] \
                                         [--profile] \
 ```                                  
for example

    python Non-Periodic/martini3-graphene-topology.py -x 21 -y 23 -z 30 -o graphene

All arguments are optional. If an argument is not used, the default value for their dimension (12 nm), and the name of the output (martini_graphene) is used.

Similarly, for generating the structure and topology of an infinite graphene sheet, run

    python Periodic/martini3-graphene-periodic.py -x [Length of the graphene sheet along x in nm] \
                                         -y [Length of the graphene sheet along y in nm] \
                                         -z [Length of the box along z in nm]
                                         -o [Name of the output for both structure and topology] \
                                         [--output-dir DIR] \
                                         [--force] \
                                         [--quiet] \
                                         [--profile] \
                                    
for example

    python Periodic/martini3-graphene-periodic.py -x 21 -y 23 -z 30 -o graphene

All arguments are optional. If an argument is not used, the default value for their dimension (12 nm), and the name of the output (martini_graphene) is used.

Both scripts also support:

- `--output-dir DIR`: write outputs to `DIR`.
- `--force`: overwrite existing output files (`.gro` and `.itp`) for the selected output name.
- `--quiet`: suppress non-essential stderr output.
- `--profile`: print stage timing information to stderr.

Without `--force`, the scripts stop if those output files already exist.

Example (overwrite existing files in a custom directory):

    python Periodic/martini3-graphene-periodic.py -x 100 -y 100 -z 30 -o graphene --output-dir results --force

## Notes
- GRO atom indices are limited to 5 digits by the file format. For systems with more than 99,999 atoms, indices wrap around.
- The two generator scripts share internal utility code from `graphene_shared.py`.

## Testing
Run the regression tests locally with:

    python -m unittest discover -s tests -p "test_*.py" -v

The same test suite is run in GitHub Actions on pushes and pull requests.

## License  
This project is licensed under the MIT License - see the [LICENSE](https://github.com/MoMS-MMSB/Martini3-Graphene/blob/main/LICENSE) file for details.

## Citation
If you use the script or the model, please cite:

Shrestha, R., Alessandri, R., Vögele, M., Hilpert, C., Souza, P. C. T., Marrink, S. J., & Monticelli, L. Martini 3 Coarse-Grained Models for Carbon Nanomaterials. *J. Chem. Theory Comput.*, 2025, 21(18), 9035–9053. https://doi.org/10.1021/acs.jctc.5c00923

<details>
<summary>BibTeX</summary>

```bibtex
@article{Shrestha2025Martini3Carbon,
  author  = {Shrestha, Roshan and Alessandri, Riccardo and V{\"o}gele, Martin and Hilpert, C{\'e}cile and Souza, Paulo C. T. and Marrink, Siewert J. and Monticelli, Luca},
  title   = {Martini 3 Coarse-Grained Models for Carbon Nanomaterials},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2025},
  volume  = {21},
  number  = {18},
  pages   = {9035--9053},
  doi     = {10.1021/acs.jctc.5c00923}
}
```

</details>
