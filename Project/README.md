README file for Multigrid Vcycle Project:

This project implements a multigrid Vcycle algorithm for solving two different partial differential equations (PDEs), namely the transient heat problem and the Poisson problem with Dirichlet boundary conditions, using Python. The main algorithm is implemented in the `Multigrid_heat.py` and `Multigrid_poison.py` files respectively.

To execute the project, run the `driver.ipynb` Jupyter Notebook file. This file imports the necessary modules and functions from the implementation files and demonstrates their use for solving the PDEs.

The `Multigrid_heat.py` file contains a class `Mgrid_transient3` which implements the multigrid Vcycle algorithm for solving the transient heat problem. The `Multigrid_poison.py` file contains a class ` Mgrid` which implements the multigrid Vcycle algorithm for solving the Poisson problem with Dirichlet boundary conditions.

The `driver.ipynb` file imports these classes and uses them to solve the respective problems. It also contains the necessary code for generating plots of the solutions.

The project is organized as follows:

```
Multigrid_vcycle_project/
├── code/
│   ├── Multigrid_heat.py
│   ├── Multigrid_poison.py
│   └── driver.ipynb
├── Images/
│   └── (contains generated plots)
└── README.md
```

The `code` folder contains the implementation files while the `data` folder is an empty directory where any data generated during the project will be saved.

To run the project, first, make sure you have all the required dependencies installed. The main dependencies for this project are:

- Python 3.x
- NumPy
- Matplotlib
- Jupyter Notebook

Once you have installed these dependencies, navigate to the `code` folder and open the `driver.ipynb` file in Jupyter Notebook. Run the cells in the notebook to execute the multigrid Vcycle algorithm for solving the transient heat and Poisson problems.

The `driver.ipynb` file contains detailed explanations of how to use the implemented classes and functions. For further details on the implementation of the multigrid Vcycle algorithm, refer to the comments in the `Multigrid_heat.py` and `Multigrid_poison.py` files.

Note: This project assumes some basic knowledge of Python, NumPy, Matplotlib, and Jupyter Notebook.

References

[1] Briggs, W. L., Henson, V. E., & McCormick, S. F. (2000). A multigrid tutorial (2nd ed.). SIAM.

[2] Zingale, M. (2017). Computational hydrodynamics for astrophysics. Retrieved from https://zingale.github.io/hydro_tutorial/

[3] Malipeddi, A. R. (2017). Geometric Multigrid Implementation. Retrieved from https://github.com/AbhilashReddyM.
