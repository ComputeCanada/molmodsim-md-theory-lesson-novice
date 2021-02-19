---
title: "Trajectory analysis with PYTRAJ"
teaching: 25
exercises: 0
questions:
- "How to analyze trajectories?"
objectives:
- "Learn how to set up and use Jupyter notebook in CC environment"
- "Learn how to visualize simulation in Jupyter notebook"
- "Learn how to calculate and plot RMSD"
- "Learn how to perform Dynamic Cross Correlation Analysis"
- "Learn how to perform Principal Component Analysis"
- "Learn how to run analysis in parallel using MPI"
keypoints:
- ""
---
### Introduction
PYTRAJ is a Python front end to the AMBER [CPPTRAJ](https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml) package. CPPTRAJ provides a variety of high level analysis commands, and at the same time it is suitable for batch processing. With PYTRAJ/CPPTRAJ you can do many operations on the raw MD trajectories. For example, convert among trajectory formats, process groups of trajectories generated with ensemble methods, image with periodic boundary conditions, create average structures, create subsets of the system. PYTRAJ is able to handle many files at the same time, and it can handle very large trajectories.

PYTRAJ offers more than 50 types of analyses such as RMS fitting, measuring distances, B-factors, radii of gyration, radial distribution functions, time correlations, and many more. PYTRAJ supports MPI, and usage of MPI is straighforward. You don't really need to understand deeply about MPI or write complicated code.

Other useful MD analysis software: [MDAnalysis](https://userguide.mdanalysis.org/stable/index.html), [Pteros](https://yesint.github.io/pteros/), [LOOS/PyLOOS](http://grossfieldlab.github.io/loos/index.htmland). These packages provide libraries that can be used to compose analysis programs. While this approach offers great flexibility, the learning curve is steep, and you will need to spend more time to master them.

References:  
1.[PTRAJ and CPPTRAJ: Software for Processing and Analysis of Molecular Dynamics Trajectory Data](https://pubs.acs.org/doi/full/10.1021/ct400341p)

### Using PYTRAJ from Jupyter Notebook

#### Installing Python Virtual Environment and Jupyter Notebook.
In this lesson we will be using PYTRAJ/AmberTools20. First you need to load modules required for AmberTools, and then load python and scipy-stack modules:

~~~
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 python scipy-stack
~~~
{: .bash}

The next step is to install and activate a virtual environment. We need virtual environment because we will be installing python modules required for this lesson.

~~~
virtualenv ~/env-pytraj
source ~/env-pytraj/bin/activate
~~~
{: .bash}

Once a virtual environment is installed and activated we can install Jupyter Notebook server. We begin installation by installing IPython kernel, the python backend for Jupyter Notebook. Kernel is a process that runs independently and interacts with the Jupyter Notebook server and its user interface. The Jupyter Notebook automatically ensures that the IPython kernel is available in the default environment. However, as we will be using Python in a specific virtualenv set up for AmberTools, we need to install IPython kernel into the newly created environment. 

~~~
pip install --no-index jupyter ipykernel
~~~
{: .bash}

To make the environment accessible from notebook we need one more step:  add the new python kernel specification to Jupyter Notebook. You can use any name fro the kernel, for example 'env-pytraj'.

~~~
python -m ipykernel install --user --name=env-pytraj
~~~
{: .bash} 

Finally, install three more packages that we will be using: 
- NGLview, a web application for molecular visualization.
- Pickle is a module providing fuctions for converting a Python object into a byte stream to store it in a file and load it back into python. 
-  Seaborn, a Python data visualization library based on matplotlib. It provides a high-level interface for drawing attractive and informative statistical graphics.

~~~
pip install seaborn pickle5 nglview
~~~
{: .bash}

For NGL viewer to work in the notebook we need to install and enable jupyter widgets extension

~~~
jupyter nbextension install widgetsnbextension --py --sys-prefix 
jupyter-nbextension enable widgetsnbextension --py --sys-prefix
~~~
{: .bash}

and enable NGLview Jupyter extension
~~~
jupyter-nbextension enable nglview --py --sys-prefix
~~~
{: .bash}

We are now ready to start Jupyter notebook server. The new Python kernel with the name `env-pytraj` will be available.


#### Launching Jupyter notebook server
This example is for launching Jupyter on Graham.

To make AmberTools available in a notebook we need to load ambertools module and activate the virtual environment before starting jupyter server. It is convenient to save all commands in a file that you can later reuse.

Let's create Jupyter startup file for use with AmberTools module: jupyter_launch_ambertools.sh with the following content: 

~~~
#!/bin/bash
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 python scipy-stack ambertools
source $EBROOTAMBERTOOLS/amber.sh
source ~/env-pytraj/bin/activate
unset XDG_RUNTIME_DIR
jupyter notebook --ip $(hostname -f) --no-browser
~~~
{: .file-content}

Before starting jupyter server we need to allocate CPUs and RAM for our notebook. We will request two MPI tasks because we will learns to how to analyze data in parallel.

~~~
salloc --mem-per-cpu=2000 --time=2:0:0 --ntasks=2
~~~
{: .bash}

~~~
salloc: Pending job allocation 44825307
salloc: job 44825307 queued and waiting for resources
salloc: job 44825307 has been allocated resources
salloc: Granted job allocation 44825307
salloc: Waiting for resource configuration
salloc: Nodes gra798 are ready for job
[svassili@gra798 ~]$ 
~~~
{:.output}

Salloc allocated the requested resources and logged you into the compute node gra798. Note the name of the node where notebook server will be running. Now we can start Jupyter server by executing commands in the file jupyter_launch_ambertools.sh

~~~
bash ./jupyter_launch_ambertools.sh
~~~
{: .bash}

Do not close this window, closing it will terminate the server. Note the port number (the default is 8888) and the notebook access token.

#### Connecting to notebook server

The message in the example above informs that notebook server is listening at port 8888 of the node gra798. Compute nodes cannot be accessed directly from the Internet, but we can connect to the login node, and the login node can connect any compute node. Thus, connection to a compute node should be also possible. How do we connect to the node gra798 at port 8888? We can instruct ssh client program to map port 8888 of gra798 to our local computer. This type of connection is called "ssh tunneling" or "ssh port forwarding". Ssh tunneling allows transporting networking data between computers over an encrypted SSH connection.

Open ANOTHER terminal tab or window and run the command:
~~~
ssh svassili@graham.computecanada.ca -L 8888:gra798:8888
~~~
{: .bash}

This SSH session created tunnel from your computer to gra798. The tunnel will be active only while the session is running. Do not close this window and do not logout, this will close the tunnel and disconnect you from the notebook.

Now in the browser you can type localhost:8888, and enter the token when prompted.

In Jupyter open new notebook. Ensure that you are creating notebook with the kernel matching the activated environment (env-pytraj), or it will fail to start!



> ## Uninstalling virtual environment from Jupyter:
>
> ~~~
> jupyter kernelspec list
> jupyter kernelspec uninstall env-pytraj
> ~~~
> {: .bash}
{: .callout}

### Computing RMSD

~~~
import pytraj as pt
import numpy as np
from matplotlib import pyplot as plt

%cd ~/scratch/Ago-RNA_sim/sim_pmemd/2-production/

# Load topology and trajectory
# iterload can load multiple trajectories, supplied as a python list
traj=pt.iterload('mdcrd',top='prmtop.parm7')

# Load reference frame
# Can also use any trajectory frame e.g ref_crd=trj[0]
ref_crd = pt.load('../../inpcrd.pdb')

# Automatically center and image molecules/residues/atoms that are outside of the box back into the box.
traj=traj.autoimage()

# Generate x-axis
tstep=0.001 # the trajectory was saved every 0.001 ns
time=np.arange(0, traj.n_frames-1)*tstep

# Compute RMSD
# Use !grep OXT ../../inpcrd.pdb to find the last residue of a protein
rmsd_data = pt.rmsd(traj, ref=ref_crd, nofit=False, mask=':1-859,@C,N,O')

# Plot rmsd
plt.plot(time,rmsd_data)
plt.xlabel("Time, ns")
plt.ylabel("RMSD, Angstrom")
~~~
{: .python}

#### Parallel trajectory analysis using MPI

~~~
%cd ~/scratch/Ago-RNA_sim/sim_pmemd/2-production/
import numpy as np
import pickle
~~~
{: .python}

Create a python script and save it as msd.py. Then submit this file to sbatch for scheduling.
~~~
%%file rmsd.py

## create a file name rmsd.py

import pytraj as pt
import pickle
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank

# load files
traj = pt.iterload(['mdcrd'], top='prmtop.parm7')
ref_crd = pt.load('../../inpcrd.pdb')

# call pmap_mpi for MPI

# we dont need to specify n_cores=6 here since we will use `mpirun -n 6`
data = pt.pmap_mpi(pt.rmsd, traj, mask=':1-859,@C,N,O', ref=ref_crd)

# data is sent to rank==0
if rank == 0:
    print(data)
    with open("rmsd.dat", "wb") as fp: 
         pickle.dump(data, fp)
~~~
{: .file-content}

Run the script on the cluster  

~~~
! srun python rmsd.py
~~~
{: .python}

Import the results pickled in the file rmsd.dat

~~~
with open("rmsd.dat", "rb") as fp: 
    rmsd=pickle.load(fp)
data=rmsd.get('RMSD_00001')
tstep=0.001
time=np.arange(0, len(data))*tstep
~~~
{: .python}

Plot RMSD

~~~
plt.plot(time,data)
plt.xlabel("Time, ns")
plt.ylabel("RMSD, Angstrom")
~~~
{: .python}


### Interactive trajectory visualization

Open a new notebook.  

- Import pytraj, nglview and make sure you are in the right directory    

~~~
import pytraj as pt
import nglview as nv
%cd ~/scratch/Ago-RNA_sim/sim_pmemd/2-production/
~~~
{: .python}   

- Load the trajectory  

~~~
traj = pt.iterload('mdcrd', top = 'prmtop.parm7')
~~~
{: .python}

- Automatically center and image molecules/residues/atoms that are outside of the box back into the box.  

~~~
traj = traj.autoimage()
~~~  
{: .python}

- Strip water and ions

~~~
trj=traj.strip(':WAT, Na+, Cl-')
~~~  
{: .python}

- Create NGLview widget

~~~
view = nv.show_pytraj(trj)
~~~  
{: .python}

- Delete the default representation

~~~
view.clear()
~~~  
{: .python}

- Add protein cartoon representation

~~~
view.add_cartoon('protein', colorScheme="residueindex", opacity=1.0)
~~~  
{: .python}

- Render view. Try interacting with the viewer using [Mouse](http://nglviewer.org/ngl/api/manual/interaction-controls.html#mouse) and [Keyboard](http://nglviewer.org/ngl/api/manual/interaction-controls.html#keyboard) controls.

~~~
view 
~~~  
{: .python}


- Add more representations. You can find samples of all representations [here](http://proteinformatics.charite.de/ngl/doc/#User_manual/Usage/Molecular_representations). 
- Try using different [coloring schemes](https://nglviewer.org/ngl/api/manual/usage/coloring.html). 
- Try visualizing different  selections. Selection language is described [here](https://nglviewer.org/ngl/api/manual/usage/selection-language.html)


~~~
view.add_licorice('protein', opacity=0.3)
view.add_hyperball(':B and not hydrogen', colorScheme="element")
view.add_hyperball(':C and not hydrogen', colorScheme="element")
view.add_spacefill('MG',colorScheme='element')
~~~  
{: .python}

- Change background color

~~~
view.background="black"
~~~  
{: .python}

- Change animation speed and step

~~~
view.player.parameters = dict(delay=0.5, step=1)
~~~  
{: .python}

- Try changing display projection

~~~
view.camera='orthographic'
~~~  
{: .python}


- Make animation smoother

~~~
view.player.interpolate = True
~~~  
{: .python}


- Set size of the widget programmatically

~~~
view._remote_call('setSize', target='Widget', args=['700px', '400px'])
~~~  
{: .python}

- Removing representations

~~~
view.remove_cartoon()
~~~  
{: .python}

- Turn on GUI

~~~
view.display(gui=True)
~~~
{: .python}

#### Useful links

- AMBER/pytraj  
  - [Atom selections](https://amber-md.github.io/pytraj/latest/atom_mask_selection.html#atom-selections)


- NGL viewer 
  - [Coloring schemes](https://nglviewer.org/ngl/api/manual/usage/coloring.html)
  - [Molecular representations](http://proteinformatics.charite.de/ngl/doc/#User_manual/Usage/Molecular_representations)
  - [Selection language](https://nglviewer.org/ngl/api/manual/usage/selection-language.html)
  - [Index](http://nglviewer.org/nglview/latest/genindex.html)
  - [Tutorial](https://ambermd.org/tutorials/analysis/tutorial_notebooks/nglview_notebook/index.html)

[Matplotlib colormaps](https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
[Seaborn colormaps](https://seaborn.pydata.org/tutorial/color_palettes.html)


#### References

1. [NGL viewer: web-based molecular graphics for large complexes](https://academic.oup.com/bioinformatics/article/34/21/3755/5021685)


### Principal component analysis
Nucleic backbone: :@O3',C3',C4',C5',O5',P


### Cross-correlation analysis. Compare windows of 200 ps. pytraj-atomiccorr
[Dynamical cross-correlation matrices](https://pubmed.ncbi.nlm.nih.gov/7563068/)

#### Jupyter Hub challenges

Bash variables can be set for a runnng python kernel from notebook:

~~~
import os
import envbash
amber_vars=os.path.join(os.environ.get("EBROOTAMBERTOOLS"), "amber.sh")
envbash.load_envbash(amber_vars)
~~~
{: .python}

I have not figured out how to load modules and set up environment variables for a python kernel running in Jupyter Hub.
