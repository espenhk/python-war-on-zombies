# Exercise E.42
import ODESolver
import numpy as np
from scitools.std import plot

class ProblemSIR:
    def __init__(self, v, beta, S0, I0, R0, T):
        '''
        v, beta: parameters in the ODE system
        S0, I0, R0: initial values
        T: simulation for t in [0, T]
        '''
        if isinstance(v, (float, int)):
            self.v = lambda t: v
        elif callable(v):
            self.v = v
            
        if isinstance(beta, (float, int)):
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta
        
        self.S0 = S0
        self.I0 = I0
        self.R0 = R0
        self.T = T
        
    def __call__(self, u, t):
        '''Right-hand side function of the ODE system.'''
        S, I, R = u
        beta = self.beta
        v = self.v
        dS = -beta(t)*S*I
        dI = beta(t)*S*I - v(t)*I
        dR = v(t)*I
        return [dS, dI, dR]

def my_beta(t):
    if t <= 12:
        return 0.0005
    else:
        return 0.0001

my_beta2 = 0.0005
my_v = 0.1
S0 = 1500
I0 = 1
R0 = 0
T = 60

my_problem = ProblemSIR(my_v, my_beta, S0, I0, R0, T)
my_problem2 = ProblemSIR(my_v, my_beta2, S0, I0, R0, T)

class SolverSIR:
    def __init__(self, problem, dt):
        self.problem, self.dt = problem, dt
        
    def solve(self, method=ODESolver.RungeKutta4):
        self.solver = method(self.problem)
        ic = [self.problem.S0, self.problem.I0, self.problem.R0]
        self.solver.set_initial_condition(ic)
        n = int(round(self.problem.T/float(self.dt)))
        t = np.linspace(0, self.problem.T, n+1)
        u, self.t = self.solver.solve(t)
        self.S, self.I, self.R = u[:,0], u[:,1], u[:,2]
        
    def plot(self):
        t = self.t
        S = self.S
        I = self.I
        R = self.R
        plot(t, S, t, I, t, R, legend=['S', 'I', 'R'])
        raw_input() # halt program to keep plot on screen, awaiting input
        
    def print_max_I(self):
        max_val = 0
        for i in self.I:
            if i > max_val:
                max_val = i
        print 'Maximum number of infected: %d' % max_val
        return max_val
        
my_dt = 0.5

solver = SolverSIR(my_problem, my_dt)
solver.solve()  # defaults to RungeKutta4
solver.plot()
max_I_1 = solver.print_max_I()    # case with beta as function printed first
solver2 = SolverSIR(my_problem2, my_dt)
solver2.solve() # see no need to plot
max_I_2 = solver2.print_max_I()   # case with beta as constant printed second
print 'Difference between first and second case: %d' % abs(max_I_1 - max_I_2)
'''
Exercise produces one plot that is with the beta as a function. If you 
would like to see the plot for beta = 0.0005, simply call solver2.plot().
Print from program follows:

Terminal> python SIR_class.py 
Maximum number of infected: 745
Maximum number of infected: 897
Difference between first and second case: 152
'''