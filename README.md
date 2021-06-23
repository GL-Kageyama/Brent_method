# Implementation comparison of Brent's method 
Brent's method is a root-finding algorithm combining the bisection method, the secant method and inverse quadratic interpolation.

## Pattern 1 : Basic implementation in Python
[0, 99.9, -19.989999999997963]  
[1, -0.04999999999999716, -4.997500000000285]  
[2, -33.33340823968873, 4434.456928842621]  
[3, -0.08746717875899111, -1.2456316167408303]  
[4, -0.09990290191258883, 0.0002707810694388968]  
[5, -0.09990019916603284, -3.360330325108407e-08]  
[6, -0.09990019950139581, 0.0]  
[-0.09990019950139581, 0.0]  

## Pattern 2 : Brent's method using Scipy's API
root: -0.09990019950139581 interation: 8  
