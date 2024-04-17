# Calligraphy Robot
This is the repo for the calligraphy robot, which implements the paper [Robot calligraphy using pseudospectral optimal control in conjunction with a novel dynamic brush model](https://arxiv.org/abs/2003.01565), published at IROS 2020. If you have any questions, please open an Issues or email sen.zephyr.wang@gmail.com

This whole program is mostly written by C++, with a few Python code to process dataset images.

## Installation dependencies
- OpenCV (both C++ and Python)
- matplotlib
- scipy
- cairosvg
- svgpathtools
- gtsam

Here are a few tips to save you some time while configuring the environment issues:
- Note, please install opencv at first. Also, if you meet problems with `imshow`, please try to solve it from the `opencv` aspect.
- Modify the `include_directories` in `CMakeLists.txt` to show the correct directory for your local `gtsam` library.


## Build and Run
To build the C++ code
```
git clone https://github.gatech.edu/borglab/calligraphy_8_1.git
cd calligraphy_8_1
mkdir build
cd build
cmake ..
make -j4
```
To generate the writing trajectory for one character, you need a copy of the SVG character [dataset](https://github.com/skishore/makemeahanzi/tree/master/svgs). Download the `svgs` folder and modify the function `separate_strokes` in `utils_node.py` to reflect that. 
```
cd ..
python3 find_poly_node.py --unicodes 2e88 --deg 3
cd build
make sen04_dynamic_optimize.run
cd ..
python showmatrix.py
```
During the optimization process, there are many parameters to adjust. To save the compile time, we read these parameters from `parameters.yaml` file. However, there are several important things that you must be careful:
- Using the same unicode for initial estimation (which is passed by `find_poly_node.py`) and optimization (which is passed by `parameters.yaml`)
- Using the same order for initial estimation (which is passed by `find_poly_node.py`) and optimization (which is passed by `Parameters.h`)
- We use absolute path in `Parameters.h`, and so you need to modify the path there to be able to read `parameters.yaml` correctly.

If you meet segmentation error, it's most probably because these parameters do not match with each other. Another possible reason is your opencv is not installed correctly that the `parameters.yaml` file is not loaded.

## Generated files explanation
All the intermediate files are stored at the `Generated_data_file` folder, for the convenience of debug.

After running the previousu code, they will generate the both the initial estimation and final optimized trajectory. The initial estimate is stored in `Generated_data_file/node_cheby_sample`. In this file, each consecutive 3 lines represent one set of Chebyshev nodes corresponding to one stroke, corresponding to x axis, y axis and z axis, respectively. x and y part use pixels as units, z part use meters. Normally, we initialize the trajectory's z component as a constant sequence, and so only the first element in the z line has a value, which is also the encoded z sequence's constant value.

The final trajectory is stored at `Generated_data_file/Sample_points_$unicode$.txt`. In this file, each line represents `x, y, z` coordinates of discrete points of the written trajectory. `x` and `y` is represented in pixel units, and `z` is represented in meters. All the points of all stroke trajectories are stored in one single file after the optimization runs.

## Results visualization
We provide many ways to visualize the final optimization results. `showmatrix.py` will show the simulation results before (`python showmatrix --file_name init`) and after (`python showmatrix --file_name opt`) optimization, and original image. `showSkeleton.py` could show the skeleton (single written trajectory). `show_z` could show the `z` trajectory after optimizing one stroke.

## Brush identification
We write a simple script `new_fit_vb_parameters.py` to do the linear regression.

# Citation
```
@article{Wang2020RobotCalligraphy,
  title={Robot Calligraphy using Pseudospectral Optimal Control in Conjunction with a Novel Dynamic Brush Model},
  author={Sen Wang and Jiaqi Chen and Xuanliang Deng and Seth A. Hutchinson and Frank Dellaert},
  journal={2020 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  year={2020},
  pages={6696-6703},
  url={https://api.semanticscholar.org/CorpusID:208157960}
}
```
