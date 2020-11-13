test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_heat2d_1 python src/heat2d.py
assert_exit_code  0
assert_in_stdout  10.        
assert_in_stdout  10.
assert_in_stdout  10.
assert_in_stdout -10.
assert_in_stdout -10.
assert_in_stdout -10.

run test_heat2d_2 python src/heat2d.py --num_elm_x 2 --num_elm_y 2 --T0_bottom 5 --T0_left -2 --heat_source 0 --flux_top 0
assert_exit_code  0
assert_in_stdout  5.        
assert_in_stdout  5.
assert_in_stdout  5.
assert_in_stdout -2.
assert_in_stdout -2.
