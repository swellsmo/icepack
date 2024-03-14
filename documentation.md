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

  - (C) Friction
    - (Cs) Side Friction
  - (u) Velocity
  - (h) Thickness
  - (τ(u,C)) Basal Shear Stress
  - (E) Energy Density
  - (D) Damage
  - (M) Membrane Stress

  ## Constants
  - year = 365.25 * 24 * 60 * 60
  - (g) gravity = 9.81 * year**2
  - (ρ_I) ice_density = 917 / year**2 * 1.0e-6
  - (ρ_W) water_density = 1024 / year**2 * 1.0e-6
  - (R) ideal_gas = 8.3144621e-3
  - (n) glen_flow_law = 3.0
  - strain_rate_min = 1e-5
  - (m) weertman_sliding_law = 3.0
  - (c) heat_capacity = 2.0e3 * year**2
  - (α) thermal_diffusivity = 2.3e-3 / (917 * 2.0) * year
  - (Tm) melting_temperature = 273.15
  - (L) latent_heat = 334e3 * year**2
  

  ## Models
  `model = icepack.models.IceShelf()`\
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
  - friction_stress(u, C):
    - Compute the shear stress for a given sliding velocity
  - bed_friction(**kwargs): 
    - Return the bed friction part of the ice stream action functional
    - velocity, friction
  - side_friction(**kwargs):
    - Return the side wall friction part of the action functional
    - velocity, thickness, side_friction
  - side_friction_xz(**kwargs): 
    - Return the side wall friction part of the action functional
    - velocity, thickness, side_friction
  - normal_flow_penalty(**kwargs):
    - Return the penalty for flow normal to the domain boundary
      - For problems where a glacier flows along some boundary, e.g. a fjord wall, the velocity has to be parallel to this boundary. Rather than enforce this boundary condition directly, we add a penalty for normal flow to the action functional.
    - velocity, scale
  
  ### Heat Transport
  Class for modeling 3D heat transport

  This class solves the 3D advection-diffusion equation for the energy
    density. The energy density factors in both the temperature and the latent
    heat included in meltwater. We use the energy density rather than the
    enthalpy because it comes out to a nice round number (about 500 MPa/m^3)
    in the unit system we use.

    surface_exchange_coefficient : float, optional
            Penalty parameter for deviation of the surface energy from the
            atmospheric value; this is very poorly constrained

Functions included:
- advective_flux(self, **kwargs):
  - energy, velocity, vertical_velocity, thickness, energy_inflow, energy_surface
- diffusive_flux(self, **kwargs):
  - energy, thickness, energy_surface
- sources(self, **kwargs):
  - energy, thickness, heat, heat_bed
- temperature(self, E):
  - Return the temperature of ice at the given energy density
- meltwater_fraction(self, E):
  - Return the melt fraction of ice at the given energy density
- energy_density(self, T, f):
  - Return the energy density for ice at the given temperature and melt fraction

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
  Transport Equation:\
  - field_name, source_name, conservative
  Continuity: Describes the form of the mass continuity equation\
  - thickness, accumulation

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
Plot a mesh with a different color for each boundary segment\
![A triplot of a sample glacier mesh, with the boundary segments labeled](https://github.com/swellsmo/icepack/assets/116534525/aa8fc611-2dd0-44ac-9847-2eafebdbf2d0)


### tricontourf
Create a filled contour plot of a finite element field\
![A filled contour plot of ice velocity using the magma colormap](https://github.com/swellsmo/icepack/assets/116534525/0ed29a91-b6a7-4b96-bf97-9a75692bb533)


### tricontour
Create a contour plot of a finite element field\
![A contour plot of the ice velocity field using the magma colormap](https://github.com/swellsmo/icepack/assets/116534525/7c0ee09e-3fc1-486b-aa47-d2634951ee41)


### tripcolor
Create a pseudo-color plot of a finite element field\
![a tripcolor plot of ice velocity](https://github.com/swellsmo/icepack/assets/116534525/08309ce2-9c8a-4920-9f79-9c234d411e59)



### quiver
Make a quiver plot of a vector field\
![A quiver plot of the model glacier's velocity using an inverted magma colormap](https://github.com/swellsmo/icepack/assets/116534525/d537ab74-db74-40b1-b314-c1bbfb6d65e9)


### streamplot
Draw streamlines of a vector field\
![A streamplot of the model glacier's velocity field in the magma colorbar](https://github.com/swellsmo/icepack/assets/116534525/fa9ef96b-1c0b-4d23-9c5b-f58946fc8d5f)



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
