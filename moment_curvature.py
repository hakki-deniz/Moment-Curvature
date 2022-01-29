"""
    This piece of code is prepared as an answer for the question at:
        https://portwooddigital.com/2022/01/16/two-fibers-five-ways/


    There are three parts of the code:
        1 - Create a function to calculate moment_curvature
        2 - Build the model and analyze it
        3 - Visualize the result
        

    The resources & references:
        1 - https://openseespydoc.readthedocs.io/en/latest/src/MomentCurvature.html
        2 - https://opensees.berkeley.edu/wiki/index.php/Moment_Curvature_Example
"""


from openseespy.opensees import *


##########
## Part 1 - Create a function
##########

def moment_curvature(integration, sec_tag, axial_load, max_d, incr_num):
    
    ## Define nodes
    node(1, 0.0, 0.0)
    node(2, 0.0, 0.0)

    ## Define boundary conditions
    fix(1, 1, 1, 1)
    fix(2, 0, 1, 0)
    
    ## Define element    
    element('zeroLengthSection',  1,   1,   2,  sec_tag)

    ## Define axial load
    timeSeries('Constant', 1)
    pattern('Plain', 1, 1)
    load(2, axial_load, 0.0, 0.0)

    ## Define analysis parameters
    integrator('LoadControl', 0.0)
    system('SparseGeneral', '-piv')
    test('NormUnbalance', 1e-4, 10) # Lesser limit than the one at the original example
    numberer('Plain')
    constraints('Plain')
    algorithm('Newton')
    analysis('Static')

    ## Create recorder
    recorder("Node", "-file", integration + ".out", "-time", "-closeOnWrite","-node", 2, "-dof", 1,2,3, "disp")
        
    ## Perform analysis -- for axial load
    analyze(1)

    ## Define reference moment to initiate analysis
    timeSeries('Linear', 2)
    pattern('Plain',2, 2)
    load(2, 0.0, 0.0, 1.0)

    ## Curvature increment
    delta_d = max_d / incr_num

    integrator('DisplacementControl', 2, 3, delta_d)

    ## Perform analysis
    analyze(incr_num)




##########
## Part 2 - Build the model
##########

model('basic','-ndm',2,'-ndf',3)

## Define material 
fy_ = 355.0      # Yield stress (MPa)
E_ = 200000.0    # Modulus of elasticity (MPa)

uniaxialMaterial('Steel01', 100, fy_, E_, 0.001) # Small strain-hardening ratio for the stress-strain relation

## Define dimensions
b_ = 300 # width of the section (mm)
d_ = 500 # depth of the section (mm)

## Variables to define fiber sections
z1 = b_/2.0
y1 = d_/2.0


#### <<< -- The exact solution -- >>>
#### The exact solution is obtained by patch command (10x10 fibers)
##integration_name = "Exact"
##section('Fiber', 1)
##patch('rect',100,10,10 ,-y1, -z1, y1, z1)

## Create fibers for a given integration
integration_name = "Radau"

section('Fiber', 1)

if integration_name == "Midpoint":
    fiber(y1/2, 0, (b_*d_/2),100)
    fiber(-y1/2, 0, (b_*d_/2),100)

elif integration_name == "Gauss":
    fiber(d_/(2*(3**(0.5))), 0, (b_*d_/2),100)
    fiber(-d_/(2*(3**(0.5))), 0, (b_*d_/2),100)

elif integration_name == "Lobatto":
    fiber(d_/2, 0, (b_*d_/2),100)
    fiber(-d_/2, 0, (b_*d_/2),100)

elif integration_name == "Radau":
    fiber(d_/2, 0, b_*d_/4,100)
    fiber(-d_/6, 0, 3*b_*d_/4,100)

elif integration_name == "Newton-Cotes":
    fiber(d_/6, 0, (b_*d_/2),100)
    fiber(-d_/6, 0, (b_*d_/2),100)


### Call the function
incr_num = 5000

target_curvature = 5**(-6)

moment_curvature(integration_name,1, 0, target_curvature, incr_num)



##########
## Part 3 - Visualize the result
##########

## Import libraries
from matplotlib import pyplot as plt
import numpy as np

## Get data from recorder output
## Exact solution
data_exact = np.loadtxt("Exact.out")  # !!! Run the analysis for the exact case first !!!

data_exact_x = data_exact[:,3]
data_exact_y = data_exact[:,0]

## Different integrations
data_ = np.loadtxt(integration_name + ".out")

data_x = data_[:,3]
data_y = data_[:,0]

## Plot the data
plt.plot(data_x, data_y, label = integration_name, linestyle = "-",color = "black")
plt.plot(data_exact_x, data_exact_y,label = "Exact", linestyle = ":",color = "black")

plt.title(integration_name + " vs Exact Solution")
plt.legend()

## Hide axis
ax = plt.gca()

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

## Show the plot
plt.show()
