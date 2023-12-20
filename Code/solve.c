// File solve.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve.h"

// Numerical derivative
double Solve_Get_Df
(double (*func)(double), double x);
// This makes the Solve functions self-contained but is not strictly necessary.

// Solve_Bisect
// Solve f(x) = nu using the bisect search method
double Solve_Bisect
(double nu, double (*func)(double), double x_min, double x_max, double tol, int *count)
{
    double x_mid, f_max, f_min, f_mid, err;
    int count_max;

    count_max = 1000; // should be large enough
    *count += 1; // adds 1 whenever Solve_Bisect is called

    // Warn and exit
    if(*count > count_max)
    {
        fprintf
        (stderr, "Solve_Bisect: Done %d iterations without convergence.\n", count_max);
        fprintf
        (stderr, "Exiting.\n");
        exit(0);
    }

    // Calculate f_max = f_(x_max) - nu
    f_max = (*func)(x_max) - nu;
    // This is how to access the function passed to this function via the function pointer (*func)

    // Calculate f_min = f_(x_min) - nu
    f_min = (*func)(x_min) - nu;

    // Warn and exit
    if(f_max * f_min > 0.0) // we cannot find a solution within the given range
    {
        fprintf
        (stderr, "Solve_Bisect: No solution within this range.\n");
        fprintf
        (stderr, "Exiting.\n");
        exit(0);
    }

    // Calculate x_mid = (x_min + x_max) / 2.0
    x_mid = (x_min + x_max) / 2.0;

    // Calculate f_mid = f_(x_mid) - nu
    f_mid = (*func)(x_mid) - nu;

    // Calculate the error
    if(nu != 0.0)
    {
        err = fabs(f_mid / nu);
    }
    else
    {
        err = fabs(f_mid);
    }

    if(err < tol) // we have a solution and the calculation ends
    {
        return x_mid;
    }

    // If f_mid * f_max < 0.0, then the solution is between x_mid and x_max
    if(f_mid * f_max < 0.0)
    {
        // Call the function Solve_Bisect with the range (x_mid, x_max)
        return Solve_Bisect(nu, func, x_mid, x_max, tol, count);
    }
    // If f_min * f_mid < 0.0, then the solution is between x_min and x_mid
    else if(f_min * f_mid < 0.0)
    {
        // Call the function Solve_Bisect with the range (x_min, x_mid)
        return Solve_Bisect(nu, func, x_min, x_mid, tol, count);
    }
    // If one of the factors is zero
    else
    {
        if (f_mid == 0.0)
        {
            return x_mid;
        }
        else if (f_max == 0.0)
        {
            return x_max;
        }
        else
        {
            return x_min;
        }
    }
}
// End of Solve_Bisect

// Solve_Get_Df
// Use f'(x) = (f(x + h) - f(x - h)) / (2h)
double Solve_Get_Df(double (*func)(double), double x_old)
{
    double h, df;

        if (x_old != 0.0)
        {
            h = x_old * 1.0E-5;
        }
        else
        {
            h = 1.0E-5;
        }

        df = (*func)(x_old + h) - (*func)(x_old - h);

        df /= 2.0 * h;
    
    return df;
}
// End of Solve_Get_Df

// Solve_Newton solves f(x) = nu using Newton's method
double Solve_Newton
(double nu, double (*func)(double), double x_0, double tol, int *count)
{
    double x_old, x_new, err, df, h;
    int count_max;

    count_max = 1000;

    // Initial value (x_0 is the starting point)
    x_old = x_0;

    do
    {
        // Do the following... while(...)

        // Get the derivative
        df = Solve_Get_Df(func, x_old);

        // Warn and exit
        if(fabs(df) < tol) // the derivative is too small
        {
            fprintf
            (stderr, "Solve_Newton: The derivative is too small.\n");
            fprintf
            (stderr, "Exiting.\n");
            exit(0);
        }
        // Warn and exit

        // Implement Newton's method
        x_new = x_old + (nu - (*func)(x_old)) / df;

        // Calculate the error
        err = fabs((x_new - x_old) / x_old);

        // Set the new x_old as x_new
        x_old = x_new;

        *count += 1;

        // Warn and exit
        if(*count == count_max) // too many iterations
        {
            fprintf
            (stderr, "Solve_Newton: Done %d iterations without convergence.\n", count_max);
            fprintf
            (stderr, "Exiting.\n");
            exit(0);
        }
    } while(err > tol); // while this condition is satisfied

    return x_new;
}
// End of Solve_Newton