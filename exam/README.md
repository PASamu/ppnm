Examination project for the course Practical Programming and Numerical Methods (Spring 2021).

Petros Adel Same
Studenternumber: 201606286

Title of the project: ODE with complex-valued functions of complex variable (Project number 20, 86 (mod 22)=20).

Description of the project: Generalize the ODE solver of your choice to solve ODE with complex-valued functions of complex variable along a straight line between the given start- and end-point.

Procedure: I use my ode-solver for real functions with real variables as my starting point. I rewrite my stepper and driver functions, such that they are suitable for differential equations with complex functions, that take complex variables. One of the biggest challenges is that instead of taking steps along a real-valued axis when doing the integration process, we now consider a whole plane because we work in the complex numbers. This problem was approached by considering absolute values. In my real ode solver I used arrays, but found it much more advantageous to work with vectors in this case.
The tools that I use to work with complex numbers are from the header files <gsl/gsl_complex.h> and <gsl/gsl_complex.math.h>.

Results: I get a segmentation fault (core dumped) when I run the program, which makes it difficult to evaluate the quality of the scripts. I am not sure where this error occurs, my best guess would be in a loop. I am quite satisfied with the rewriting of the step and driver functions such that they work for complex differential equations. When I look through the code, I take the same steps and go through the same equations as for the real case. 

Final rating of project: 5/10.
