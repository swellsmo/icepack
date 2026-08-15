
# icepack [![DOI](https://zenodo.org/badge/45697304.svg)](https://zenodo.org/badge/latestdoi/45697304)

icepack is a library for modeling the flow of ice sheets and glaciers using the finite element method.
For more information and installation instructions, see the [project webpage](https://icepack.github.io).
Once you've installed icepack, the directory `notebooks/tutorials/` contains several Jupyter notebooks with demonstrations of the main features of icepack to get you started.
The directory `notebooks/how-to` contains more advanced material for when you've got a good feel for how things work.
Some of these demos use real data, for which you'll need to have an account with [NASA EarthData](https://urs.earthdata.nasa.gov/).

See also [firedrake](https://www.firedrakeproject.org), a finite element modeling package used to implement the solvers in icepack.


# Documentation for non-computery people
Maintained by Sarah Wells-Moran (swellsmo@uchicago.edu)

Disclaimer: This is my best interpretation of how icepack works. It could be very, very wrong lol

A lot of this is probably unnecessary or over explained. My purpose in making and maintaining this documentation is to create the document I wish I had when I was starting out as an undergrad and had 0 idea how computers worked. Even as a grad student, I have maybe 3 brain cells that understand what computering is. 

Hardware/software in case it's relevant to anyone: 2024 Mac Mini with an M4 chip running Tahoe 26.6.1<sup>[[1]](#tahoe)</sup>

## Docker 🐳

Docker is a containerization software. Containerization is the concept of sending out software packages with all of the files, code, and libraries it needs to run straight out of the box on any operating system. Below, I will be talking about Docker images and containers. A docker image is a read-only template that tells the daemon how to create containers with all of the dependencies needed to run things installed. The container is a running instance of the image. You can have multiple containers running off of a single image, but the containers do not talk to each other or the files in the host machine. 

---

  <details open>  
    <summary> <b>Dockerfile</b> </summary>
    The Dockerfile is the cornerstone of the docker process. It tells the Docker daemon (the daemon manages everything behind the scenes) which components to include in your image build. 
 
  ```
# check=skip=InvalidBaseImagePlatform
FROM firedrakeproject/firedrake-vanilla:2025-01

ENV PORT=8870

RUN sudo apt update && sudo apt install patchelf
RUN . /home/firedrake/firedrake/bin/activate && \
    git clone https://github.com/icepack/icepack.git && \
    python3.12 -m pip install --upgrade pip && \
    pip3 install --editable ./icepack && \
    pip install ipykernel && \
    pip3 install --upgrade jupyter jupyterlab notebook jupyter_server jupyter_core && \
    pip3 install netcdf4 ipympl && \
    python3 -m pip install siphash24 && \
    python3 -m ipykernel install --user --name=firedrake

EXPOSE "$PORT"

# Set up a shell script within the container to start the virtual environment 
# and open jupyter notebooks when entering container
USER 0

COPY entry.sh /usr/local/bin/
RUN echo -e "\njupyter notebook --ip 0.0.0.0 --no-browser --port $PORT" >> /usr/local/bin/entry.sh && \
    chmod +x /usr/local/bin/entry.sh

# Set up bashrc for interactive sessions (courtesy of Justin P Linick)
RUN touch /root/.bashrc && \
    echo "PS1='🐳 \[\e[1;32m\]\u@ice-fem\[\e[m\]:\[\e[1;34m\]\w\[\e[m\]\\$ '" >> /root/.bashrc

USER firedrake

ENTRYPOINT ["/usr/local/bin/entry.sh"]
```

<details> 
  <summary> Breaking down the Dockerfile line by line </summary>
  
* `# check=skip=InvalidBaseImagePlatform` is a Parser Directive and is optional. This specific directive suppresses an (in my opinion) annoying warning that pops up if you build the image on an operating system that is different from the operating system the base image was built from. I've included more details below for anyone who is interested
     * <details> <summary> What is a parser directive? </summary>
       When running a build, Docker passes the Dockerfile through a parser, which then hands the contents of the file over to the daemon. A Parser Directive tells the parser how to handle the contents of the file. A common example of this is setting the escape character, which tells the daemon to treat the character as a new line. On Linux and Unix platforms, the escape character is `\`, but in Windows, that denotes a file path. You can add a parser directive `# escape=\` to tell the parser to treat that character as the escape character. </details>
     * <details> 
       <summary> What does the parser directive 'check' do? </summary> 
       Docker runs build checks when building an image to see if the Dockerfile adheres to their predefined "Best Practices". Following the rules helps avoid errors, but also can be kind of annoying if you've found ways to mitigate the errors elsewhere in your workflow. The `# check` parser directive tells the parser how to run build checks on the image. In this case, the check skips the step of checking which build platform (operating system) the base image uses versus which build platform your machine is running. A full list of build checks is available [here](https://docs.docker.com/reference/build-checks/). When a check fails, it doesn't mean the build has failed, but it might throw up errors or fail in future versions of Docker.  </details>
     * Note: If including this, it MUST go in the first line of the Dockerfile, even before any comments. Otherwise, the parser will interpret it as another comment.
* `FROM firedrakeproject/firedrake-vanilla:2025-01` is the base image you want to build your Docker image on top of. More details about the `FROM` command are available [here](https://dockerbuild.com/reference/from)
* `ENV PORT=8870` Sets the environment variable PORT to 8870. Environmental variables are available to containers running off of this image. 
     * I use this variable so I can change the port my computer communicates to the container on just in one location. I specifically use ENV here over ARG because my start.sh script accesses the PORT variable to determine which port to publish in the `docker run` command.
* `RUN` commands run the following bits of code. Each command is separated by `&&`, and `\` denotes a line break to the parser. The only two worth mentioning here that aren't package installation are:
     * `. /home/firedrake/firedrake/bin/activate`: This activates the firedrake virtual environment and must be the first command of the `RUN` block. From my limited understanding, firedrake is installed inside the virtual environment, not outside of it. Thus, we must first activate the virtual environment and install everything icepack needs inside of the environment. We will also need to activate the virtual environment any time we want to use icepack (more on this later). Further reading recommended on the icepack install page: <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>
     * `python3 -m ipykernel install --user --name=firedrake`: This makes a kernel called 'firedrake' that Jupyter notebooks can communicate with when we run the demo notebooks. This must be run inside the virtual environment
* `EXPOSE "$PORT"`: Tells Docker that the system listens on the specified network port. Here, I'm using the predefined variable `$PORT` so that I only need to change the value in one location if I want to expose a different port.
     * NOTE: This does not publish the port. Publishing the port means setting up a connection between the Docker container and the computer. EXPOSE is mostly a communication to whoever is running the container that the specific port needs to be published. [This](https://dhavalgojiya.hashnode.dev/understanding-dockers-expose-keyword-4-port-mapping-scenarios-explained) write up provides some good examples of what `EXPOSE` means in terms of container behavior.
* `USER 0`: Sets the user to root so that we can mess around with things in the root folder without running into permissions errors
* `COPY entry.sh /usr/local/bin/` This takes a file entry.sh from my build context and copies its contents to /usr/local/bin within the image. entry.sh is a small shell script I built so I don't have to repeatedly type `source ~/firedrake/bin/activate` when opening the container lol. I have included my build context folder in the github, but the shell script is relatively simple:
     * <details> 
       <summary> entry.sh</summary>

        ```
        #!/bin/bash

        source ~/firedrake/bin/activate
        cd ~/icepack/notebooks
        ```
        </details>

* `RUN echo -e "\njupyter notebook --ip 0.0.0.0 --no-browser --port $PORT" >> /usr/local/bin/entry.sh && \`
     * This adds a new line to the entry.sh file that reads `jupyter notebook --ip 0.0.0.0 --no-browser --port 8870`
          * The --port flag is necessary for the notebook to run and communicate with your computer! 
     * I add this line here rather than in the script file itself because I want to be able to change the port value via the `$PORT` variable
* `chmod +x /usr/local/bin/entry.sh`: continuation of the previous `RUN` command. Changes permissions on the entry.sh script so that any user can run it, not just the root user
* OPTIONAL: `RUN touch /root/.bashrc && \`: Creates a file called .bashrc in the /root directory.
* OPTIONAL:  `echo "PS1='🐳 \[\e[1;32m\]\u@firedrake\[\e[m\]:\[\e[1;34m\]\w\[\e[m\]\\$ '" >> /root/.bashrc`
     * Adds the line `PS1='🐳 \[\e[1;32m\]\u@firedrake\[\e[m\]:\[\e[1;34m\]\w\[\e[m\]\\$ '` to /root/.bashrc.
     * .bashrc is a file that the Bash shell reads and runs every time you start up an interactive session. More info available [here](https://www.digitalocean.com/community/tutorials/bashrc-file-in-linux)
     * The variable `PS1` controls the prompt within your Bash terminal session. This particular one turns your username green and adds a little whale emoji beside it :)
     * I don't typically see this behavior because I'm in the notebook and my entry.sh script bypasses command line stuff
* `USER firedrake`: changes user from 0 (root) back to firedrake. It is best practice to limit root access in applications, and we don't really need it anyways because we're not changing root directory aside from
* `ENTRYPOINT ["/usr/local/bin/entry.sh"]`: Tells Docker to run the script entry.sh (that we copied and modified earlier in the script) every time a container starts
     * The script entry.sh activates the firedrake virtual environment, changes our directory to the ~/icepack/notebooks folder, and starts up a jupyter notebook. Basically we enter the container and boom we immediately get a link to open the notebooks. Change the `cd ~/icepack/notebooks` line if the notebooks you're working on are in a different directory
     * Entrypoint can be overridden in the `docker run` command with the flag `--entrypoint <command>`. When I am troubleshooting the image, I often run with `--entrypoint /bin/bash` so I can poke around inside the virtual environment without having the notebook up.
     * Also worth it to note, under this setup, the container will end as soon as you quit jupyter. It works best this way for my workflow but everyone has different preferences! You can modify entry.sh within the build context folder to change what happens at startup.

</details> <!--- End of breaking down dockerfile --->

</details> <!--- End of Dockerfile section --->

---

  <details open>
    <summary> <b>Building a Docker image</b> </summary>
Create a file named Dockerfile in your machine and paste the above text into it. Then in terminal, navigate to the folder that contains Dockerfile and run:
    
```
  docker build -t <image name>:<tag> .
```
    
#### Tag
The `<tag>` feature allows you to have multiple instances of a Docker image with the same name. For example, if you have an image named `icepack` and you want to modify it slightly but keep the old instance in case things break, you can use the command `docker build -t icepack:main .` for your base image and `docker build -t icepack:1.0.0 .` for your modified image. It is important to keep old versions of your images so you can revert back to something you know works when you inevitably break things! 

#### Build Context
The `.` at the end of the command is the Build Context. This tells the builder what files it has access to when building. In this case, `.` tells Docker to build the image with all the folders and files in the current directory. If the dockerfile is in a different directory, replace `.` with the relative path to the dockerfile directory. For example, I keep my Dockerfile in a subdirectory called 'docker' within my icepack folder because otherwise it tries to access the .git directory in my icepack folder and throws up errors. When I build, I set the context as `./docker`. I found this website to be a super good overview of how the build context works: <https://dockerbuild.com/reference/build-context>.

#### Other flags
`-f, --file`: By default, the docker build command looks for the file `PATH/Dockerfile` to build the image off of. If your Dockerfile has any extensions or name other than `Dockerfile`, e.g. Dockerfile.txt or Dockerfile_icepack, the builder won't be able to find the file. Add the flag `-f` to the command to direct it to the location and name of the file you want to build off of:
```
docker build -t <image name>:<tag> -f /path/to/Dockerfile.txt .
```
You can run this command from any directory in terminal if you use the absolute path to the Dockerfile. The Dockerfile does not need to be in the build context directory when you use this flag. However, running this command from any directory can cause a lot of lag in building your Docker image as it will send the entire build context to the daemon to build the image. 


#### More information
I stumbled across [this website](https://dockerbuild.com/courses/) while trying to find the best words to explain the build context in an easy-to-understand way. It does a phenomenal job at breaking everything down even further than the Docker manuals, which is great if you're like me and totally new to Docker. If you're interested in learning more about how a Dockerfile works and how to build your own, the Courses tab is a great place to start! 

  </details> <!--- End of Building a Docker Image --->

> [!NOTE]
> <a name="tahoe"><b>[1]</b></a>: For some reason upgrading to Tahoe breaks a lot of things about Docker (I wish I had stuck with Sequoia). In particular, Tahoe breaks websocket connections with Docker, which is the open line of communication between your Docker container and the Jupyter Notebook. I don't know what changed, but I started getting stuck in websocket_timeout loops when trying to do anything on my notebooks. Thankfully, there was a simple fix: go to Docker Desktop, open settings and look for `General>Choose Virtual Machine Manager (VMM)` and make sure the toggle for Use Rosetta is OFF. 

  ## Starting Shell Script
  I use a shell script to open the docker image and immediately start the jupyter notebook. This allows me to run a simple command rather than having to remember the exact flags and mounts I need within the container. This script is based off of my workflow and file organization, so you'll need to modify the mounts to your liking. I've included a detailed description of what everything does below. 

  
  Create a file <name of script>.sh. Paste the code block below into the sh file. Open the terminal and run `sh <name of script>.sh`. 
  ```
  #!/bin/bash
set -e

IMAGE_NAME="<image name>:<tag>" # Change this to your project name (must be lowercase)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIR_ENV_DIR="/home/firedrake"

mkdir -p "$SCRIPT_DIR/notebooks"
mkdir -p "$SCRIPT_DIR/notebooks/data"
mkdir -p "$SCRIPT_DIR/notebooks/meshes"
mkdir -p "$SCRIPT_DIR/src"

MOUNTS=(
    -v "$SCRIPT_DIR:$VIR_ENV_DIR/icepack"
    -v "$SCRIPT_DIR/notebooks:$VIR_ENV_DIR/icepack/notebooks"
    -v "$SCRIPT_DIR/notebooks/data:$VIR_ENV_DIR/icepack/notebooks/data"
    -v "$SCRIPT_DIR/notebooks/meshes:$VIR_ENV_DIR/icepack/notebooks/meshes"
    -v "$SCRIPT_DIR/src:$VIR_ENV_DIR/icepack/src"
)

docker run --rm -it \
    --platform linux/amd64 \
    "${MOUNTS[@]}" \
    -p 8887:8887 \
    "$IMAGE_NAME"
```
### Breaking down the components of the shell script:
`set -e` stops the execution of a script if any command in the script has an error and returns the exit code of the command within the script that failed

`IMAGE_NAME`, `SCRIPT_DIR`, and `VIR_ENV_DIR` are variables used by the script. 

> [!TIP]
> At least on Mac, you can't open the Jupyter notebooks in browser unless you have the `-p` or `--publish` flag in your docker run command. The above script publishes the docker container to port 8887, which is the same port exposed in the Dockerfile

  ## Vocaublary
  *args - non-keyworded arguments (Requires correct order of arguments passed to a function)\
  $\textsf{\color{#3F00FF}**kwargs}$ - Keyword arguments (Arguments can be passed in any order, in the style of kwarg=argument)

  - (C) Friction
    - (Cs) Side Friction
  - (u) Velocity
  - (h) Thickness
  - (τ(u, C)) Basal Shear Stress
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
  The model tells the solvers what physics problem we want to solve. The solvers then input the relevant arguments into the model and tell us the solution for the physics problem given the paramters specified. 
  
  Pass $\textsf{\color{#3F00FF}**kwargs}$ into the solver for each model, not the model itself. All of the $\textsf{\color{#3F00FF}**kwargs}$ presented in this section are the $\textsf{\color{#3F00FF}**kwargs}$ that can be input to the solvers of each model

  The full list of models is:
  ```DamageTransport
  HeatTransport3D
  HybridModel
  IceShelf
  IceStream
  ShallowIce
  ```
  
### Damage Transport
```
damage_model = icepack.models.DamageTransport()
damage_solver = icepack.solvers.DamageSolver(damage_model)
```
  Description of the continuum damage mechanics model
  
  This module contains a solver for the conservative advection equation that describes the evolution of ice damage (Albrecht and Levermann 2014).
  

  This model contains the following predefined values (that can be changed with $\textsf{\color{#3F00FF}**kwargs}$):
  - damage_stress [0.07]
  - damage_rate [0.3]
  - healing_strain_rate [2e-10 * year]
  - healing_rate [0.1]

  Note: The source code does not allow for the above kwargs to be added, I modified the code to do this

  Additional $\textsf{\color{#3F00FF}**kwargs}$: damage, strain_rate, membrane_stress, damage_inflow, flux,
  
  
### Heat Transport
```
heat_model = icepack.models.HeatTransport3D()
heat_solver = HeatTransportSolver(heat_model)
```
Class for modeling 3D heat transport

This class solves the 3D advection-diffusion equation for the energy density. The energy density factors in both the temperature and the latent heat included in meltwater. We use the energy density rather than the enthalpy because it comes out to a nice round number (about 500 MPa/m^3) in the unit system we use.

    surface_exchange_coefficient : float, optional
            Penalty parameter for deviation of the surface energy from the
            atmospheric value; this is very poorly constrained

The functions defined in this model are used by the Heat Transport Solver. All the $\textsf{\color{#3F00FF}**kwargs}$ listed here can be changed within the solver.

Functions defined in model:
- advective_flux(self, **kwargs):
  - $\textsf{\color{#3F00FF}**kwargs}$: energy, velocity, vertical_velocity, thickness, energy_inflow, energy_surface
- diffusive_flux(self, **kwargs):
  - $\textsf{\color{#3F00FF}**kwargs}$: energy, thickness, energy_surface
- sources(self, **kwargs):
  - $\textsf{\color{#3F00FF}**kwargs}$: energy, thickness, heat, heat_bed
- temperature(self, E):
  - Return the temperature of ice at the given energy density
- meltwater_fraction(self, E):
  - Return the melt fraction of ice at the given energy density
- energy_density(self, T, f):
  - Return the energy density for ice at the given temperature and melt fraction

### Hybrid
Can resolve both plug and shear flow
```
model = icepack.models.HybridModel(friction=friction_function)
solver = icepack.solvers.FlowSolver(model, **opts)
```

Functions defined outside of model:
- gravity(**kwargs)
  -  velocity, thickness, surface
- terminus(**kwargs)
  - Return the terminal stress part of the hybrid model action functional
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness, surface
  - $\textsf{\color{red}**opts}$: ice_front_ids
-  _effective_strain_rate(ε_x, ε_z, ε_min)
- stresses(**kwargs)
  - Calculate the membrane and vertical shear stresses for the given horizontal and shear strain rates and fluidity
  - $\textsf{\color{#3F00FF}**kwargs}$: strain_rate_x, strain_rate_z, fluidity, strain_rate_min
- horizontal_strain_rate(**kwargs)
  - Calculate the horizontal strain rate with corrections for terrain-following coordinates
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, surface, thickness
- vertical_strain_rate(**kwargs)
  - Calculate the vertical strain rate with corrections for terrain-following coordinates
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness
- viscosity(**kwargs)
  - Return the viscous part of the hybrid model action functional
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, surface, thickness, fluidity, strain_rate_min

$\textsf{\color{#3F00FF}**kwargs}$ taken by model:
- viscosity
- friction
- side_friction
- side_friction_xz
- gravity
- terminus
- penalty
- continuity

$\textsf{\color{red}**opts}$ taken by model: 
- ice_front_ids
- side_wall_ids

Functions defined inside model:
- action(self, **kwargs)
  - Return the action functional that gives the hybrid model as its Euler-Lagrange equations
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness
  - $\textsf{\color{red}**opts}$: ice_front_ids, side_wall_ids
- scale(self,**kwargs)
  - Return the positive, convex part of the action functional. The positive part of the action functional is used as a dimensional scale to determine when to terminate an optimization algorithm.
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity
- quadrature_degree(self,**kwargs)
  - Return the quadrature degree necessary to integrate the action functional accurately
    - Firedrake uses a very conservative algorithm for estimating the number of quadrature points necessary to integrate a given expression. By exploiting known structure of the problem, we can reduce the number of quadrature points while preserving accuracy.
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness

### Ice Shelf
Class for modelling the flow of floating ice shelves
```
model = icepack.models.IceShelf()
solver = icepack.solvers.FlowSolver(model, **opts)
```

This class provides functions that solve for the velocity and thickness of a floating ice shelf. The relevant physics can be found in ch. 6 of Greve and Blatter.
  
Functions defined outside of the model:
- gravity(**kwargs)
  - Return the gravitational part of the ice shelf action functional
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness
- terminus(**kwargs)
  - Return the terminus stress part of the ice shelf action functional
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness


$\textsf{\color{#3F00FF}**kwargs}$ taken by model:
- viscosity
- gravity
- terminus
- side_friction
- penalty
- continuity

$\textsf{\color{red}**opts}$ taken by model:
- ice_front_ids
- side_wall_ids

Functions defined inside of model:
- action(self, **kwargs)
  - Return the action functional that gives the ice shelf diagnostic model as the Euler-Lagrange equations
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness
  - $\textsf{\color{red}**opts}$: ice_front_ids, side_wall_ids
- scale(self, **kwargs)
  - Return the positive, convex part of the action functional
  - The positive part of the action functional is used as a dimensional scale to determine when to terminate an optimization algorithm.
- quadrature_degree(self, **kwargs)
  - Return the quadrature degree necessary to integrate the action functional accurately
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness

  
### Ice Stream
Assumes plug flow, e.g. the velocity is roughtly constant with depth
```
model = icepack.models.IceStream(friction=friction_function)
solver = icepack.solvers.FlowSolver(model, **opts)
```

Functions defined outside of model:
- gravity(**kwargs)
  - Return the gravitational part of the ice stream action functional: $E(u) = -\int_\Omega\rho_Igh\nabla s\cdot u\; dx$
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness, surface
- terminus(**kwargs)
  - Return the terminal stress part of the ice stream action functional
  - The power exerted due to stress at the ice calving terminus $\Gamma$ is $E(u) = \frac{1}{2}\int_\Gamma\left(\rho_Igh^2 - \rho_Wgd^2\right) u\cdot \nu\, ds$
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness, surface


$\textsf{\color{#3F00FF}**kwargs}$:
- viscosity
- friction
- side_friction
- penalty
- gravity
- terminus
- continuity

### Shallow Ice
Assumes shear flow, e.g. the speed at the ice base is much smaller than at the ice surface

## Other items in the "Models" folder:
These documents are in the models folder of icepack, but do not contain any models to solve. Instead, they contain functions that other models call. Here is a list of defined functions in each document, as well as the $\textsf{\color{#3F00FF}**kwargs}$ they can take, which can be input into the Solver

### Friction

  Functions included:
  - friction_stress(u, C):
    - Compute the shear stress for a given sliding velocity and friction
  - bed_friction(**kwargs): 
    - Return the bed friction part of the ice stream action functional
    - $\textsf{\color{#3F00FF}**kwargs}$: velocity, friction
  - side_friction(**kwargs):
    - Return the side wall friction part of the action functional
    - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness, side_friction
  - side_friction_xz(**kwargs): 
    - Return the side wall friction part of the action functional
    - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness, side_friction
  - normal_flow_penalty(**kwargs):
    - Return the penalty for flow normal to the domain boundary
      - For problems where a glacier flows along some boundary, e.g. a fjord wall, the velocity has to be parallel to this boundary. Rather than enforce this boundary condition directly, we add a penalty for normal flow to the action functional.
    - $\textsf{\color{#3F00FF}**kwargs}$: velocity, scale

### Transport
Transport Equation:
- field_name, source_name, conservative
  
Continuity: Describes the form of the mass continuity equation
- thickness, accumulation

### Viscosity
Constants included:
- transition_temperature = 263.15  [K]
- A0_cold = 3.985e-13 * year * 1.0e18  [1 / (MPa^3 yr)]
- A0_warm = 1.916e3 * year * 1.0e18  [1 / (MPa^3 yr)]
- Q_cold = 60  [kJ / mol]
- Q_warm = 139  [kJ / mol]

Functions included: 
- rate_factor(T)
  - Compute the rate factor in Glen's flow law for a given temperature
- _effective_strain_rate(ε, ε_min)
- membrane_stress(**kwargs)
  - Calculate the membrane stress for a given strain rate and fluidity
  - $\textsf{\color{#3F00FF}**kwargs}$: strain_rate, fluidity, strain_rate_min
- viscosity_depth_averaged(**kwargs)
  - Return the viscous part of the action for depth-averaged models
  - $\textsf{\color{#3F00FF}**kwargs}$: velocity, thickness, fluidity, strain_rate_min

## Solvers


IcePack contains 5 solvers

```
FlowSolver
ImplicitEuler
LaxWendroff
HeatTransportSolver
DamageSolver
```

The main solvers I use are FlowSolver, HeatTransportSolver, and DamageSolver



### Damage Solver
Solver for the continuum damage mechanics model, damage_transport()
  
### Flow Solver
When initiating a flow solver, pass it any arguments ($\textsf{\color{red}**opts}$) that never change throughout a simulation.\
`solver = icepack.solvers.FlowSolver(model, **opts)`

Includes ImplicitEuler and LaxWendroff solvers

##### $\textsf{\color{red}**opts}$ can include:
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

#### Diagnostic Solver

#### Prognostic Solver
Incluldes a timestep component, dt, that must be the first position in the arguments\
accumulation might be involved somehow?

### Heat Transport Solver
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






List of color variables:\
blue=#3F00FF\
function=#7F00FF\
