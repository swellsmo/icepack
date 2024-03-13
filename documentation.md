# Documentation
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
  **kwargs - Keyword arguments (Arguments can be passed in any order)

  ## Models
  `model = icepack.models.IceShelf()`
  
  ### Damage Transport
  *Unsure if I can pass **kwargs in this one, need to test
  
  variable_name [default value]\
  damage_stress [0.07]\
  damage_rate [0.3]\
  healing_strain_rate [2e-10 * year]\
  healing_rate [0.1]

  ### Friction
  
  ### Heat Transport

  ### Hybrid

  ### Ice Shelf
  Class for modelling the flow of floating ice shelves

  This class provides functions that solve for the velocity and
  thickness of a floating ice shelf. The relevant physics can be found
  in ch. 6 of Greve and Blatter.
  
  velocity\
  thickness\
  viscosity\
  gravity\
  terminus\
  side_friction\
  penalty
  
  ### Ice Stream

  ### Shallow Ice

  ### Transport
## Flow Solver
  When initiating a flow solver, pass it any arguments that never change throughout a simulation.\
  `solver = icepack.solvers.FlowSolver(model, **opts)`\
  **opts can include things like dirichlet_ids, side_wall_ids, ice_front_ids

  ### Diagnostic Solver

  ### Prognostic SOlver









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
