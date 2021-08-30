# Opening and Closing Surfaces
Public code release for ["Opening and Closing Surfaces"](http://dgp.toronto.edu/~sgsellan/pdf/opening-and-closing-surfaces.pdf), presented at SIGGRAPH Asia 2020 and authored by [Silvia Sellán](http://dgp.toronto.edu/~sgsellan/), Ang Yan Sheng, Jacob Kesten and [Alec Jacobson](http://www.cs.toronto.edu/~jacobson/). Released under MIT License.

## Installation
We have only succesfully validated this installation on MacOS 10.15 with Matlab 2020a, but are confident it should work on most Unix-based machines and versions of Matlab from the past five years. 

To install this library, please start by cloning the repository recursively
```
git clone --recursive https://github.com/sgsellan/opening-and-closing-surfaces.git
```
After this, we will build the mex functions in the `gptoolbox` directory:
```
cd opening-and-closing-surfaces/gptoolbox/mex
mkdir build
cd build
cmake .. -DMatlab_ROOT_DIR=“MATLAB_ROOT_DIRECTORY” (i.e. /usr/local/MATLAB/R2021a)
make
```
Then, build our own mex files by entering Matlab and, in Matlab, adding `matlab-include` recursively to your Matlab path (for instance, by running `addpath(genpath(path/to/opening-and-closing-surfaces/matlab-include))`), navigating to `matlab-include/mex/` and running `build_mex` in the Matlab console.

## Use
We allow you to replicate results from our paper exactly, by runnning the scripts in the `figures/` directory. Instructions to run and understand the output of each can be found in each scripts' first commented lines. In most cases, the scripts will generate input and output `.obj` files which we then rendered using Blender for the paper figures. You can look at our Blender setup in `render/` and substitute the existing meshes with the input and output from our scripts to exactly replicate our paper figures up to very minor lighting direction and orientation choices.

## Known Issues
Please do not hesitate to contact
[sgsellan@cs.toronto.edu](mailto:sgsellan@cs.toronto.edu) if you find any issues
or bugs in this code, or you struggle to run it in any way.

## Graphics Replicability Stamp Initiative
We, the authors, hereby give the Graphics Replicability Stamp Initiative and its reviewers permission to review the code and advertise the review publicly after the stamp is approved.


