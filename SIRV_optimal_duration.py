# Exercise E.45
'''
Program more or less equivalent to E.44, except from where explicitly noted.
'''

import ODESolver
import numpy as np
from scitools.std import plot

class ProblemSIR:
    def __init__(self, v, beta, p, S0, I0, R0, V0, T):
        '''
        v, beta, p: parameters in the ODE system
        S0, I0, R0, V0: initial values
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
        
        if isinstance(p, (float, int)):
            self.p = lambda t: p
        elif callable(p):
            self.p = p
        
        self.S0 = S0
        self.I0 = I0
        self.R0 = R0
        self.V0 = V0
        self.T = T
        
    def __call__(self, u, t):
        '''Right-hand side function of the ODE system.'''
        S, I, R, V = u
        beta = self.beta
        v = self.v
        p = self.p
        dS = -beta(t)*S*I - p(t)*S
        dI = beta(t)*S*I - v(t)*I
        dR = v(t)*I
        dV = p(t)*S
        return [dS, dI, dR, dV]

class SolverSIR:
    def __init__(self, problem, dt):
        self.problem, self.dt = problem, dt
        
    def solve(self, method=ODESolver.RungeKutta4):
        self.solver = method(self.problem)
        ic = [self.problem.S0, self.problem.I0, self.problem.R0, 
              self.problem.V0]
        self.solver.set_initial_condition(ic)
        n = int(round(self.problem.T/float(self.dt)))
        t = np.linspace(0, self.problem.T, n+1)
        u, self.t = self.solver.solve(t)
        self.S, self.I, self.R, self.V = u[:,0], u[:,1], u[:,2], u[:,3]
        
    def plot(self):
        t = self.t
        S = self.S
        I = self.I
        R = self.R
        V = self.V
        plot(t, S, t, I, t, R, t, V, legend=['S', 'I', 'R', 'V'])
        raw_input() # halt program to keep plot on screen, awaiting input
        
    def print_max_I(self):
        max_val = 0
        for i in self.I:
            if i > max_val:
                max_val = i
        print 'Maximum number of infected: %d' % max_val
        return max_val

t_array = np.linspace(0, 31, 31+1)
I_array = np.zeros(31+1)
for Vt in range(0, 31+1):
    print 'Vt=%d' % Vt,
    my_beta = 0.0005
    my_v = 0.1
    def my_p(t):
        if 6 <= t <= 6 + Vt:
            return 0.1
        else:
            return 0
    S0 = 1500
    I0 = 1
    R0 = 0
    V0 = 0  # assume 0 vaccinated at start of scenario
    T = 60

    my_problem = ProblemSIR(my_v, my_beta, my_p, S0, I0, R0, V0, T)
        
    my_dt = 0.5

    solver = SolverSIR(my_problem, my_dt)
    solver.solve()  # defaults to RungeKutta4
    max_I = solver.print_max_I()
    I_array[Vt] = max_I
plot(t_array, I_array)

'''
Terminal> python SIRV_optimal_duration.py 
Vt=0 Maximum number of infected: 877
Vt=1 Maximum number of infected: 764
Vt=2 Maximum number of infected: 669
Vt=3 Maximum number of infected: 591
Vt=4 Maximum number of infected: 532
Vt=5 Maximum number of infected: 490
Vt=6 Maximum number of infected: 462
Vt=7 Maximum number of infected: 447
Vt=8 Maximum number of infected: 441
Vt=9 Maximum number of infected: 441
Vt=10 Maximum number of infected: 441
Vt=11 Maximum number of infected: 441
Vt=12 Maximum number of infected: 441
Vt=13 Maximum number of infected: 441
Vt=14 Maximum number of infected: 441
Vt=15 Maximum number of infected: 441
Vt=16 Maximum number of infected: 441
Vt=17 Maximum number of infected: 441
Vt=18 Maximum number of infected: 441
Vt=19 Maximum number of infected: 441
Vt=20 Maximum number of infected: 441
Vt=21 Maximum number of infected: 441
Vt=22 Maximum number of infected: 441
Vt=23 Maximum number of infected: 441
Vt=24 Maximum number of infected: 441
Vt=25 Maximum number of infected: 441
Vt=26 Maximum number of infected: 441
Vt=27 Maximum number of infected: 441
Vt=28 Maximum number of infected: 441
Vt=29 Maximum number of infected: 441
Vt=30 Maximum number of infected: 441
Vt=31 Maximum number of infected: 441
'''
'''
The plot shows us that the graph approaches a constant (~441), but our
prints reveal that the maximum number of infected is constant when
Vt >= 8. We can also see that the decrease from Vt=7 to Vt=8 is quite
small, but considering the context of human lives it seems worth it
to vaccinate for one extra day to save these six lives.
'''