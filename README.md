ThinPlateSplines.jl
--------

This small module implements [thin-plate splines](https://en.wikipedia.org/wiki/Thin_plate_spline) for Julia. Thin-plate splines are a very useful technique for representing smooth, continuous deformations in 2, 3, or N-dimensional space. They are often used in data interpolation, shape matching, and geometric design. Thin-plate splines allow you to specify a deformation on the entire space by specifying a set of 'control points' - their original and deformed positions. You can then 'extrapolate' the resulting deformation to the rest of the space, in a way that minimizes the bending energy over the entire space.

In this module, deformations are represented using the `ThinPlateSpline` type. Three functions are exposed:

1. `tps_solve`, which takes as input a set of control points and deformed control points, and a stiffness coefficient, and produces a deformation. 
2. `tps_deform`, which takes as input as set of points and a deformation, and deforms the points according to the deformation.
3. `tps_energy`, a utility function which returns the bending energy of a thin-plate spline. This function is useful when you want to compare different deformations to each other.

#### Example usage
Consider the following points (vertices of a triangle in 2D):

```julia
x1 = [0.0 1.0 
      1.0 0.0
      1.0 1.0]
```

And consider deforming them slightly:

```julia
x2 = [0.0 1.0
      1.1 0.0
      1.2 1.5]
```

Now calculate the TPS associated with this deformation (setting the stiffness coefficient to 1.0):

```julia
tps = tps_solve(x1,x2,1.0)
```

And deform another set of points using this deformation (note that one of the points is the same, and will be deformed in the same way):

```julia
x = [1.0 0.0
     2.0 2.0]
tps_deform(x,tps)
```

The output should be something like:
```
 1.1  5.55112e-17
 2.5  3.5
```

We can also calculate the bending energy of this deformation:

```
tps_energy(tps)
```

The result should be 0, which indicates that it is a fully affine transform.

#### Transforming lines/curves

How about a more general example? How about deforming *continuous lines*, not just points? `tps_deform` is general enough to be able to do that, without having to actually discretize the line as a lot of points. We can represent lines parametrically using linear functions of some variable `t`, with the [TaylorSeries](https://github.com/JuliaDiff/TaylorSeries.jl) package. For example, let's represent the line in the plane given by `x(t) = 2 + t, y(t) = 1 + t/2`:
```julia
using TaylorSeries
t = Taylor1([0,1],15)
line = [2 + t, 1 + 0.5*t]'
```

We can deform this line just as we can with points:
```julia
deformed_line = tps_deform(line, tps)
```

This will look something like: [`2.4 + 1.25 t + O(t¹⁶)`,  `2.0 + 1.25 t + O(t¹⁶)`]. In this case, since the transform is fully affine, the result is also a straight line, and so no terms higher than `t` exist. In general tps transforms, higher-order terms may exist. For instance:

```julia
x1 = [0.0 1.0
      1.0 0.0
      1.0 1.0
      1.0 2.0]
x2 = [0.0 1.0
      1.1 0.1
      1.2 1.5
      1.3 2.0]
tps_deform(line,tps_solve(x1,x2,1.0))
```

Gives a more complicated curve, with `y(t) = 1.472 + 0.603 t - 0.009 t² + ...`e