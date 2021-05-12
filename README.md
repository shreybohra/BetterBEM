# BetterBEM
A fast Blade Element Method algorithm with guaranteed convergence, based off the equations reduced by Ning (2014)

## How To Use

The dcp_calc function outputs the contribution to Cp and forces experienced by an element of the blade.
The input parameters are as follows:

1. **twist** - the twist of the blade element
2. **polars** - a cell array containing the polars of the airfoil used at the element - see below for more details
3. **z** - the radial position of the element
4. **R** - the maximum radius of the blade
5. **tsr** - the tip speed ratio of the blade
6. **B** - the total number of blades used in the completed design
7. **chord** - the chord of the aerofoil used 
8. **r0** - the distance from the centre of rotation to the root of the blade
9. **V** - incoming fluid velocity

## Aerofoil Polars
Polars is a cell array containing four different fits.
Elements 1 and 2 are fits of the lift coefficient and drag coefficient respectively, against angle of attack **in radians**, between -10 and 30 degrees. This can be obtained from XFOIL.

Elements 3 and 4 are fits of Cl and Cd respectively, but between angles of attack betwene -180 degrees and +180 degrees. 
