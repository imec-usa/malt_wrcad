***  Configuration File  ***
*circuit file name: netlist.cir
*parameter file name: 
*passfail file name: 
*envelope file name: 

***  Circuit Node Defaults  ***
dt = 4e-11
dx = 2

***  Circuit Parameter Defaults  ***
nominal = 1
min = 0.5
max = 2
sigma = 0
sig_abs = 0
static = 0
logs = 1
corners = 0
include = 1

***  Circuit Nodes  ***
node = v(phi.Xb0), dt = 4e-11, dx = 2
node = v(phi.Xb0.XI1), dt = 4e-11, dx = 2
node = v(phi.Xb1.XI12), dt = 4e-11, dx = 2

***  Circuit Parameters  ***
param = Xlcomp, nominal = 1, min = 0.3, max = 3, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 1, include = 0
param = Xjcomp, nominal = 1, min = 0.3, max = 3, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 1, include = 0
param = Xac, nominal = 1, min = 0.3, max = 1.7, sigma = 4, sig_abs = 0, static = 0, logs = 0, corners = 0, include = 1
param = Xpdc, nominal = 1, min = 0.3, max = 1.7, sigma = 4, sig_abs = 0, static = 0, logs = 0, corners = 0, include = 1, nom_min = 1, nom_max = 1
param = inoff, nominal = 0, min = -1, max = 1, sigma = 0, sig_abs = 0.04, static = 0, logs = 0, corners = 0, include = 1
param = inipp, nominal = 1, min = 0, max = 3, sigma = 0, sig_abs = 0.04, static = 0, logs = 0, corners = 0, include = 1
param = b0, nominal = 0.08, min = 0.03, max = 0.3, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 0, include = 1
param = l0, nominal = 5, min = 1.5, max = 15, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 0, include = 1
param = b1, nominal = 0.018, min = 0.006, max = 0.06, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 0, include = 0
param = b2, nominal = 0.07, min = 0.02, max = 0.2, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 0, include = 0
param = b34, nominal = 0.04, min = 0.012, max = 0.12, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 0, include = 0
param = l45, nominal = 30, min = 10, max = 100, sigma = 4, sig_abs = 0, static = 0, logs = 1, corners = 0, include = 0

***  2D Margins  ***
***  Example. Replace arguments with actual param names  ***
* param_x = XJ, param_y = XL

***  File Extensions  ***
circuit_extension = .cir
parameters_extension = .param
passfail_extension = .passf
envelope_extension = .envelope
plot_extension = .plot

***  General Options  ***
*fraction of sigma:
binsearch_accuracy = 0.1
spice_call_name = wrspice
max_subprocesses = 0
threads = 16
verbose = 0
print_terminal = 1

***  Options for Define  ***
d_simulate = 1
d_envelope = 1

***  Options for Margins  ***

***  Options for 2D Margins  ***
2D_iter = 16

***  Options for Corners Yield  ***
***  ranges: 0-10, 0-9, 1-40
y_search_depth = 5
y_search_width  = 5
y_search_steps  = 12
y_max_mem_k = 4194304
y_accuracy = 0.1
y_print_every = 0

***  Options for Optimize  ***
o_min_iter = 100
o_max_mem_k = 4194304
