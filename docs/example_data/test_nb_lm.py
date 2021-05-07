import numpy as np
import nlopt
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')



# Tutorial problem
import nlopt
from numpy import *

def myfunc(x, grad):
    if grad.size > 0:
        grad[0] = 0.0
        grad[1] = 0.5 / sqrt(x[1])
    return sqrt(x[1])

def myconstraint(x, grad, a, b):
    if grad.size > 0:
        grad[0] = 3 * a * (a*x[0] + b)**2
        grad[1] = -1.0
    return (a*x[0] + b)**3 - x[1]

opt = nlopt.opt(nlopt.LN_COBYLA, 2)
opt.set_lower_bounds([-float('inf'), 0])
opt.set_min_objective(myfunc)
opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,2,0), 1e-8)
opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,-1,1), 1e-8)
opt.set_xtol_rel(1e-4)
x = opt.optimize([1.234, 5.678])
minf = opt.last_optimum_value()
print("optimum at ", x[0], x[1])
print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())


# My problem with NLppt
z = np.random.poisson(np.linspace(50, 100, 100))
z = np.array([185.11374566, 202.45933565, 204.80998553, 203.83234955,
       200.21618182, 212.18184362, 185.33562681, 216.69946935,
       227.77438551, 234.26695672, 215.85697455, 151.25482012,
       162.29698644, 152.03128138, 121.93018198, 110.22190348,
       119.45182617,  12.82944133,  50.30388017,   4.53817323])
z = np.array([ 27.05387981,  40.65447387, 126.21675658, 121.19161715,
       109.50943745,  31.20712814,  47.82935946, 101.73260216,
       110.27419303,  65.74332385, 115.33171395,  77.01233731,
        94.0765563 ,  80.55924159,  95.74064971,  68.67075162,
         5.29780232,   1.76400602,   1.77175652,   7.06877443])
plt.scatter(range(z.size), z, s=10)

mu = 1
def func(x, grad):
    size = z.size
    k_vec = np.arange(size) + 1
    res = 0
    for i in range(size):
        lam = x[0] + (i+1)*x[1]
        res = res + z[i]*(np.log(lam) - np.log(lam + mu)) - mu*np.log(lam + mu)
    return res

opt = nlopt.opt(nlopt.LN_BOBYQA, 2)
opt.set_lower_bounds([1e-3, -float('inf')])
opt.set_max_objective(func)
opt.set_xtol_rel(1e-4)
x = opt.optimize([np.mean(z), 0.])

minf = opt.last_optimum_value()
print("optimum at ", x)
print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())

?opt

np.linspace(1, 10)
observations = np.random.poisson(np.linspace(10, 1, 10))
k_vec = np.arange(observations.size) + 1
mu = 1
def func(params):
    alpha, beta = params
    res = 0
    for i in range(k_vec.size):
        lam = alpha + beta*(i+1)
        res = res + observations[i]*(np.log(lam + 1e-6) - np.log(lam + mu+ 1e-6)) - mu*np.log(lam + mu+ 1e-6)
    return -res

initial_guess = [1, -1]
result = optimize.minimize(func, initial_guess)
if result.success:
    fitted_params = result.x
    print(fitted_params)
else:
    raise ValueError(result.message)

result

def func2(observations, mu, lam):
    size = observations.size
    k_vec = np.arange(size) + 1
    res = 0
    for i in range(size):
        res = res + observations[i]*(np.log(lam) - np.log(lam + mu)) - mu*np.log(lam + mu)
    return res

lam = np.linspace(5, 15, 100)

y = func2(observations, 1, lam)
plt.plot(lam, y)

alphas = np.linspace(5, 15, 20)
betas = np.linspace(0, 2, 2)

# Should have ML estimates alpha=10, beta=0
observations = np.random.poisson(10, size=100)


def f(params):
    # print(params)  # <-- you'll see that params is a NumPy array
    a, b, c = params # <-- for readability you may wish to assign names to the component variables
    return a**2 + b**2 + c**2

initial_guess = [1, 1, 1]
result = optimize.minimize(f, initial_guess)
if result.success:
    fitted_params = result.x
    print(fitted_params)
else:
    raise ValueError(result.message)


np.random.poisson(np.linspace(100))
plt.scatter(range(observations.size), observations)
observations
func2(observations, 1, 5)

alphas.shape
betas.shape
X,Y = np.meshgrid(alphas, betas) # grid of point
Z = func(observations, 1, X, Y) # evaluation of the function on the grid
im = plt.pcolor(Z) # drawing the function
# adding the Contour lines with labels
plt.colorbar(im) # adding the colobar on the right
plt.xlabel('alpha')
plt.ylabel('beta')
# latex fashion title
plt.show()

np.arange(10)
