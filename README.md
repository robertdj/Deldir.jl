# Deldir

[![Build Status](https://travis-ci.org/robertdj/Deldir.jl.svg?branch=master)](https://travis-ci.org/robertdj/Deldir.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/ox6gslc6nm58sbka?svg=true)](https://ci.appveyor.com/project/robertdj/deldir-jl)
[![codecov.io](https://codecov.io/github/robertdj/Deldir.jl/coverage.svg?branch=master)](https://codecov.io/github/robertdj/Deldir.jl?branch=master)

[deldir](https://cran.r-project.org/package=deldir) is an R package for computing Delaunay triangulations and Voronoi/Dirichlet tesselations.
This Julia package wraps the Fortran code from deldir.


## Usage

The coordinates of the generators are specified as two vectors that are fed to `deldir`, the main function of this package:
```julia
x = rand(8)
y = rand(8)
del, vor, summ = deldir(x, y)
```

The output from `deldir` are three [DataFrames](https://github.com/JuliaData/DataFrames.jl):
One for the topology of the Delaunay triangulation; one for topology of the Voronoi tesselation; one with a summary mainly related to the area of the triangles and Voronoi cells.

By default, `deldir` works with points in the unit rectangle, but other bounding rectangles can be specified as a third argument.

The area of the Voronoi cells are also available directly with the function `voronoiarea`.

Two functions are available to extract the edges of the Delaunay triangles and Voronoi cells in a "plot friendly" manner:
```julia
Dx, Dy = edges(del)
Vx, Vy = edges(vor)
```

Using the results from above this can be plotted using e.g. the [Plots package](https://github.com/tbreloff/Plots.jl):

```julia
using Plots
scatter(x, y, xlim = (0,1), ylim = (0,1), markersize = 6, label = "generators")
plot!(Dx, Dy, label = "Delaunay")
plot!(Vx, Vy, style = :dash, label = "Voronoi")
```

One realization looks like the following.

![Delaunay & Voronoi edges](deldir.png)


## Installation

Install the package by running

```julia
]add Deldir
```


# Details

There are other Julia package for interacting with Voronoi cells and Delaunay triangulations.
The deldir R package has been around for a long time and therefore it is my hope that the Fortran code give correct results.
Furthermore, deldir has two qualities I value:

- It interacts well with the bounding box.
- It returns the area of the Voronoi cells *in the same order as the input generators*.

I am also the author of the pure Julia package [VoronoiCells](https://github.com/JuliaGeometry/VoronoiCells.jl) with similar functionalities.
The *VoronoiCells* package executes *much* faster.
Consider the time taken to run the `voronoiarea` functions of both packages with an increasing number of points:

![Comparison of Deldir and VoronoiCells](comparison.png)

The script generating this output is available in the `examples` folder.
The comparison plot is made with
```julia
julia> versioninfo()
Julia Version 1.1.0
Commit 80516ca202* (2019-01-21 21:24 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, skylake)
```


## Compiled code

To make it easier to install *Deldir* the Fortran code is cross-compiled using the [BinaryBuilder package](https://github.com/JuliaPackaging/BinaryBuilder.jl) in a [dedicated repository](https://github.com/robertdj/DeldirBuilder).

I am not using the latest version of the Fortran code, because the cross-compilation fails on macOS from version 0.1-16 and onwards.
I do not have access to a contemporary Mac, so it is difficult for me to troubleshoot the issue.
If you are using a Mac and would like newer features you are welcome to submit a pull request.
The [changelog](https://cran.r-project.org/web/packages/deldir/ChangeLog) for the R package may provide guidance to what happened.


## Limitations

Not all features of the R package are available.
I have e.g. chosen to ignore options regarding dummy points. 

Are you missing anything important? 
Check out the [manual](https://cran.r-project.org/web/packages/deldir/deldir.pdf) for the R package to see if the Fortran library supports it.


# Acknowledgement

[Rolf Turner](https://www.stat.auckland.ac.nz/~rolf) is author of the deldir package for R as well as all Fortran code in this package.


# License

The *Julia code* in this package is MIT licensed and the *Fortran code* is licensed under GPL.

