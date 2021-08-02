# python-war-on-zombies
A basic epidemiology analysis assignment to calculate the spread of a viral pandemic, including vaccination policies and zombies waking from the dead.
 
Coded in the 1st semester of college, as part of the INF1100 (today IN1900) course at the University of Oslo.

These are all either book exercises or mandatory assignments for that course. 
The book exercises can be found in "A Primer on Scientific Programming with Python,
3rd Edition" by H.P. Langtangen. Other versions of the book will likely also work,
but I can not guarantee the exercise numbers will match up.

## Files:
* `ODESolver.py` is a framework for solving (Ordinary Differential Equations), with
  manually written implementations of popular ODE solving algorithms like 
  RungeKutta4, Forward Euler and Backward Euler. It also imports a method that 
  implements Newton's Method for ODE solving. 
* `SIR.py` and `SIR_class.py`: Exercises E.35 and E.36, where the latter is making
  the former solution more flexible by creating a class to solve the problem.
* `SIRV.py`: Exercise E.37: builds on `SIR.py`, by adding vaccination.
* `SIRV_varying_p.py`: Exercise E.38: alternative to `SIRV.py`, by implementing a vaccination
  campaign.
* `SIRV_optimal_duration.py`: Exercise E.39: builds on `SIRV_varying_p.py` by 
  looking for the optimal vaccination intensity, start point and duration.
* `SIZR.py`: Exercise E.40: introduces **Z**ombies to the SIR problem, where infected
  that die go on to become zombies.
* `Night_of_the_Living_Dead.py`: Exercise E.41, modelling the entire movie "Night
  of the Living Dead", by extending the SIZR problem.
* `test.py` has test functions for most of these assignments.