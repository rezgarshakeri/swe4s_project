# Solving 2D Heat Equation with Finite Element Method

This repository contains a Finite Elements solver for 2D Heat Equation.

## How to run
This code can be run with:

`python src/heat2d.py`

If you want to change the number of mesh in x and y directions, or change the value of boundary conditions you can run with:

`python heat2d.py --num_elm_x 4 --num_elm_y 4 --T0_bottom 0 --T0_left -5 --heat_source 5 --flux_top 10`

num_elm_x: number of element in x direction

num_elm_y: number of element in y direction

T0_bottom: prescribed temperature on the bottom edge

T0_left: prescribed temperature on the left edge

heat_source: Magnitude of the source in the domain

flux_top: heat flux on the top boundary

This returns the solution array which represents the Temperature of all nodes in the discrete domain. 