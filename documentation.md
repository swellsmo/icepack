# Documentation
Disclaimer: This is my best interpretation of how icepack works. It could be very, very wrong lol
  ## Commands to Open Docker Container
  Build a container:\
  `docker build -t icepack_image:0.1.4 -f F:\icepack\Dockerfile.txt .`

  Run the container:\
  `docker run -it --rm --publish 8888:8888 -v F:\icepack:/products -v "F:\QGIS Layers":/mynotebooks/data -v F:\icepack\Notebooks:/mynotebooks icepack_image:0.1.3`\
  `source ~/firedrake/bin/activate`\
 `cd /mynotebooks`\
  `jupyter notebook --ip 0.0.0.0 --no-browser`

  ## Vocaublary
  *args - non-keyworded arguments (Requires correct order of arguments passed to a function)\
  **kwargs - Keyword arguments (Arguments can be passed in any order, in the style of kwarg=argument)

  ## Models
  `model = icepack.models.IceShelf()`
  All models need velocity, thickness
  
  ### Damage Transport
  *Unsure if I can pass **kwargs in this one, need to test?
  
  variable_name [default value]\
  damage_stress [0.07]\
  damage_rate [0.3]\
  healing_strain_rate [2e-10 * year]\
  healing_rate [0.1]\

  **kwargs\
  damage\
  strain_rate\
  membrane_stress\
  damage_inflow

  ### Friction
  Functions included:
  - bed_friction
  - side_friction
  - side_friction_xz
  - normal_flow_penalty
  
  ### Heat Transport

  ### Hybrid
  Can resolve both plug and shear flow

  surface\

  Functions Included:
  - _effective_strain_rate(ε_x, ε_z, ε_min)
  - stresses(**kwargs)
  - horizontal_strain_rate(**kwargs)
  - vertical_strain_rate(**kwargs)

  ### Ice Shelf
  Class for modelling the flow of floating ice shelves

  This class provides functions that solve for the velocity and
  thickness of a floating ice shelf. The relevant physics can be found
  in ch. 6 of Greve and Blatter.
  
  viscosity\
  gravity\
  terminus\
  side_friction\
  penalty
  
  ### Ice Stream
  Assumes plug flow, e.g. the velocity is roughtly constant with depth
  

  ### Shallow Ice
  Assumes shear flow, e.g. the speed at the ice base is much smaller than at the ice surface

  ### Transport

  ### Viscosity
  Functions included: 
  - rate_factor(A)
  - _effective_strain_rate(ε, ε_min)
  - membrane_stress(**kwargs)
  - viscosity_depth_averaged(**kwargs)
    -  velocity, thickness, fluidity, strain_rate_min

  ## Damage Solver
  Solver for the continuum damage mechanics model

  flux\
  damage
  
  ## Flow Solver
  When initiating a flow solver, pass it any arguments that never change throughout a simulation.\
  `solver = icepack.solvers.FlowSolver(model, **opts)`
  #### **opts can include:
        model: 
            The flow model object -- IceShelf, IceStream, etc.
        dirichlet_ids : list of int, optional
            Numerical IDs of the boundary segments where the ice velocity
            should be fixed
        side_wall_ids : list of int, optional
            Numerical IDs of the boundary segments where the ice velocity
            should have no normal flow
        diagnostic_solver_type : {'icepack', 'petsc'}, optional
            Use hand-written optimization solver ('icepack') or PETSc SNES
            ('petsc'), defaults to 'icepack'
        diagnostic_solver_parameters : dict, optional
            Options for the diagnostic solver; defaults to a Newton line
            search method with direct factorization of the Hessian using
            MUMPS
        prognostic_solver_type : {'lax-wendroff', 'implicit-euler'}, optional
            Timestepping scheme to use for prognostic equations, defaults
            to Lax-Wendroff
        prognostic_solver_parameters : dict, optional
            Options for prognostic solve routine; defaults to direct
            factorization of the flux matrix using MUMPS

  ### Diagnostic Solver

  ### Prognostic Solver
  Incluldes a timestep component, dt, that must be the first position in the arguments\
  accumulation might be involved somehow?

  ## Heat Transport Solver
  aflux - advective flux\
  dflux - diffusive flux\
  sources\
  energy\
  thickness\



## Graph Types

### triplot
Plot a mesh with a different color for each boundary segment

### tricontourf
Create a filled contour plot of a finite element field

### tricontour
Create a contour plot of a finite element field

### tripcolor
Create a pseudo-color plot of a finite element field

### quiver
Make a quiver plot of a vector field

### streamplot
Draw streamlines of a vector field


## Creating a Mesh


  # Dockerfile
  ```
  FROM firedrakeproject/firedrake-vanilla:2023-11
CMD ["/bin/bash"]
RUN sudo apt update && sudo apt install patchelf
RUN . /home/firedrake/firedrake/bin/activate && \
    pip install git+https://github.com/icepack/Trilinos.git && \
    pip install git+https://github.com/icepack/pyrol.git && \
    git clone https://github.com/swellsmo/icepack.git && \
    pip install ./icepack && \
    python -m ipykernel install --user --name=firedrake && \
    pip install --upgrade jupyter jupyterlab jupyter_server jupyter_client jupyter_core && \
    pip install netcdf4
EXPOSE 8888
ENV VIRTUAL_ENV ~/firedrake
ENV PATH "$VIRTUAL_ENV/bin:$PATH"
CMD ["/bin/bash"]
```
