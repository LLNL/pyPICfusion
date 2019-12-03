# pyPICfusion

pyPICfusion is a set of python scripts to perform nuclear fusion for an arbitrary distribution of ions. The scripts were developed to be used to demonstrate and validate algorithms as a part of particle-in-cell (PIC) code framework. This work has been published: "A Pairwise Nuclear Fusion Algorithm for Weighted Particle-in-Cell Plasma Simulationsi" D. P. Higginson, A. Link, A. Schmidti. Journal of Computational Physics **388**, 439 (2019)[doi:10.1016/j.jcp.2019.03.020](https://doi.org/10.1016/j.jcp.2019.03.020).

## Running the code

The scripts are divided into three parts.

* **PARTMAKER**: initializes the ion distributions that will be used for fusion
* **FUSION**: performs the fusion interactions between the ions creates the fusion products
* **MAKEPLOTS**: performs different analysis routines on the fusion products

Each of these parts has a "inputdeck.py" and a "source.py". In general, the source portion contains functions that are relatively unchanged. The inputdecks will be changed by the user to modify the parameters of interest. 

### PARKMAKER

This portion of the code initialized the ion distributions. The output of this is an hdf file that contains all of the ion reactant particles desired. The user may want to create their own ion particles, for instance from the output of another code. Then the use can create these particles and pass an object to PARTMAKER\_inputdeck.save\_particle\_reactants(S).

Here "S" is an object with the following attributes, with N being the number of ion particles

* S.W  -- a (N,) numpy array of the the particle weights
* S.vx -- a (N,) numpy array of the velocities in the x-direction
* S.vy -- a (N,) numpy array of the velocities in the y-direction
* S.vz -- a (N,) numpy array of the velocities in the z-direction

### FUSION

This script controls the fusion of the ions. It requires the ions generated in the previous step.

### MAKEPLOTS

This script performs analysis of the output fusion products. 

## Authors

* **Drew P. Higginson**, Lawrence Livermore National Laboratory

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

LLNL-CODE-759242

