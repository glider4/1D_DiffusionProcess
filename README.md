# 1D_DiffusionProcess

Estimating a 1D heat equation diffusion process via Explicit, Implicit, and Crank-Nicolson methods.  These 
solutions have been implemented in NumPy using Linear Algebra.  The results themselves are hard to tell apart. 
The Heat Equation itself can be modelled in 1D, 2D, or 3D to show the diffusion of heat across a number of different 
constructs.  Here, it is only estimated in 1D for general use and to show the Linear Algebra solution.

One advantage of using Linear Algebra here is the loop's overall size in terms of code is fairly short - for each of these
estimation schemes, < 10 lines would be needed to run the algebraic formuals through as long as it is set up properly 
outside of the loop.  Compare this to some of my other scripts which have lengthy estimation properties, these are a 
significant decrease in coding necessary.  However, one drawback is the readability of `np.diag` and `np.eye`, etc. can 
be more convoluted for generally programming use outside of the world of applied mathematicians.

Regardless, the results show for themselves the estimation power that these schemes have.  These plots were the result of 
running this script with `n = 20000` and `dx = 0.01`.

![result](./plots.png)
