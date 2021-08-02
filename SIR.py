# Exercise E.41
from numpy import linspace
from ODESolver import RungeKutta4
from scitools.std import plot, hold

'''
S + I + R = N
S' + I' + R' = 0
S(t + dt) = S(t) - beta*S*I*dt
S'(t) = -beta*S*I
I(t + dt) = I(t) + beta*S*I*dt - v*I*dt
I'(t) = beta*S*I - v*I
R(t + dt) = R(t) + v*I*dt
R'(t) = v*I
'''
S0 = 1500
I0 = 1
R0 = 0

def my_terminate(u, t, step_no):
    '''
    Terminate function to be used with the ODESolver environment.
    Takes in a n x 3 - array and checks if the current values of
    S, I and R for the given step add up (within a tolerance) to
    the same as the original S0, I0 and R0 put together. If this
    tolerance is breached, the function returns True and stops the
    solution.
    '''
    S = u[:, 0]
    I = u[:, 1]
    R = u[:, 2]
    
    if abs((S[step_no] + I[step_no] + R[step_no]) - (S0 + I0 + R0)) > 1:
        return True
    else:
        return False
        
def solve_problem(beta, v):
    '''
    Function to solve the ODE problem in E.41. Takes in a beta and
    v value, sets some initial variables and calculates the number
    of steps to the solution, creates a t_array for the time points
    on which to solve the ODE. Then defines a function f to return
    the system of ODEs in the exercise, gives this function to an
    instance of the RungeKutta4 class, sets initial conditions and
    solves the ODE on the t_array with the my_terminate function as
    a condition to keep solving.
    Returns the two arrays u_sol and t_sol, where u_sol is a 
    n x 3 - array.
    '''
    dt = 0.5
    t_min = 0
    t_max = 60
    steps = int(round(float(t_max - t_min)/dt))
    t_array = linspace(t_min, t_max, steps)
    U0 = [S0, I0, R0]
    
    def f(u, t):        
        S, I, R = u         # u given as 3-tuple
        dS = -beta*S*I
        dI = beta*S*I - v*I
        dR = v*I
        return [dS, dI, dR]
    
    solver = RungeKutta4(f)
    solver.set_initial_condition(U0)
    u_sol, t_sol = solver.solve(t_array, my_terminate)
    return u_sol, t_sol
    
def plot_solution(u, t):
    '''
    Function takes in a n x 3 - array from the solve_problem-function,
    assigns the variable names S, I and R to the respective columns of
    the array and plots these against the t array. Halts the program,
    awaiting keyboard input to remove the window (and f.ex. proceed to
    the next plot).
    '''
    S = u[:, 0]
    I = u[:, 1]
    R = u[:, 2]
    plot(t, S, t, I, t, R, legend=['S', 'I', 'R'])
    raw_input()
    
# Variables to be used in our two cases:
beta1 = 0.0005
beta2 = 0.0001
v = 0.1

# calls to functions, solving the exercise:
u_sol1, t_sol1 = solve_problem(beta1, v)
plot_solution(u_sol1, t_sol1)
u_sol2, t_sol2 = solve_problem(beta2, v)
plot_solution(u_sol2, t_sol2)

'''
We can clearly see that where in the first case of beta = 0.0005 nearly
all susceptibles have become infected (or even recovered) by day 20 (and
definitely by day 30), staying inside (beta = 0.0001) means less than 100
people are removed from the susceptibles group in the first 60 days.
This is definitely a drastic improvement, which teaches us all to stay inside
during epidemics.
'''