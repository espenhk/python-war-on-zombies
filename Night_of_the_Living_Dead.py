# Exercise E.47
import ODESolver
import numpy as np
from scitools.std import plot

class ProblemSIR:
    def __init__(self, sigma, beta, d_S, d_I, p, alpha, S0, I0, Z0, R0, T):
        '''
        sigma, beta: parameters in the ODE system
        S0, I0, Z0, R0: initial values
        T: simulation for t in [0, T]
        '''
        if isinstance(sigma, (float, int)):
            self.sigma = lambda t: sigma
        elif callable(sigma):
            self.sigma = sigma
            
        if isinstance(beta, (float, int)):
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta
            
        if isinstance(d_S, (float, int)):
            self.d_S = lambda t: d_S
        elif callable(d_S):
            self.d_S = d_S
            
        if isinstance(d_I, (float, int)):
            self.d_I = lambda t: d_I
        elif callable(d_I):
            self.d_I = d_I
            
        if isinstance(p, (float, int)):
            self.p = lambda t: p
        elif callable(p):
            self.p = p
            
        if isinstance(alpha, (float, int)):
            self.alpha = lambda t: alpha
        elif callable(alpha):
            self.alpha = alpha
        
        self.S0 = S0
        self.I0 = I0
        self.Z0 = Z0    # exercise specific
        self.R0 = R0
        self.T = T
        
    def __call__(self, u, t):
        '''Right-hand side function of the ODE system.'''
        S, I, Z, R = u  # Z added
        beta = self.beta
        d_S = self.d_S
        d_I = self.d_I
        p = self.p
        alpha = self.alpha
        sigma = self.sigma
        # dS, dI, dZ and dR formulas changed
        dS = sigma(t) - beta(t)*S*Z - d_S(t)*S
        dI = beta(t)*S*Z - p(t)*I - d_I(t)*I
        dZ = p(t)*I - alpha(t)*S*Z
        dR = d_S(t)*S + d_I(t)*I + alpha(t)*S*Z
        return [dS, dI, dZ, dR] # dZ added to return

class SolverSIR:
    def __init__(self, problem, dt):
        self.problem, self.dt = problem, dt
        
    def solve(self, method=ODESolver.RungeKutta4):
        self.solver = method(self.problem)
        ic = [self.problem.S0, self.problem.I0, self.problem.Z0,
              self.problem.R0]  # Z0 added
        self.solver.set_initial_condition(ic)
        n = int(round(self.problem.T/float(self.dt)))
        t = np.linspace(0, self.problem.T, n+1)
        u, self.t = self.solver.solve(t)
        # self.Z added:
        self.S, self.I, self.Z, self.R = u[:,0], u[:,1], u[:,2], u[:,3]
        
    def plot(self):
        t = self.t
        S = self.S
        I = self.I
        Z = self.Z
        R = self.R
        plot(t, S, t, I, t, Z, t, R, legend=['S', 'I', 'Z', 'R']) # plot of Z
        raw_input()
        
    def print_max_I(self):
        max_val = 0
        for i in self.I:
            if i > max_val:
                max_val = i
        print 'Maximum number of infected: %d' % max_val
        return max_val

# Following code is exercise-specific:
def my_beta(t):
    if t <= 4:
        return 0.03
    elif 4 < t <= 28:   # first 4 hrs + 24 hrs = 28
        return 0.0012
    else:
        return 0
        
def my_alpha(t):
    if t <= 4:
        return 0
    elif 4 < t <= 28:
        return 0.0016
    else:
        return 0.006
        
def my_d_I(t):
    if t <= 4:
        return 0
    elif 4 < t <= 28:   # first 4 hrs + 24 hrs = 28
        return 0.0014
    else:
        return 0

def my_d_S(t):
    if t <= 4:
        return 0
    elif 4 < t <= 28:
        return 0
    else:
        return 0.0067
        
def my_sigma(t):
    if t <= 4:
        return 20
    elif 4 < t <= 28:
        return 2
    else:
        return 0
        
my_p = 1 # constant for all phases, so no need to define function

# initial conditions, and as such not functions of time:
my_S0 = 60
my_Z0 = 1
my_I0 = 0
my_R0 = 0
my_T = 4 + 24 + 5

my_dt = 0.5

my_zombie_problem = ProblemSIR(my_sigma, my_beta, my_d_S, my_d_I, 
                               my_p, my_alpha, my_S0, my_I0, my_Z0, 
                               my_R0, my_T)
my_zombie_solution = SolverSIR(my_zombie_problem, my_dt)
my_zombie_solution.solve()
my_zombie_solution.plot()

'''
Exercise produces one plot that is the SIZR development. 
Nothing is printed.
Note that population appears to increase - this is expected behaviour, as
our sigma is greater than 0 and constitutes the amount of people brought
into the zombified area.
'''