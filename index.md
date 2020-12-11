
## benchmark monitoring 


### intersection benchmarks

The following benchmarks test the surface intersections, the x-axis are commits in time, so want to see a negative gradient.

#### This tests shows the intersections with all sensitive surfaces of the TrackML detector:

It performs 10000 tracks with each track trying to intersect all surfaces int he TrackML detector without preselection.

![TML Intersections](figures/TML_INTERSECT_ALL.png)
![All Surfaces Benchmark](figures/BM_INTERSECT_ALL.png)

#### This tests shows the intersection with concentric cylinders
Cylinders are positioned at (0,0,0) and not rotated wrt the global frame
![Concentric Cylinder Benchmark](figures/BM_INTERSECT_CONCETRIC_CYLINDERS.png)

#### This tests shows the intersection with arbitrary cylinders
![Generic Cylinder Benchmark](figures/BM_INTERSECT_CYLINDERS.png)

#### This tests shows the intersection with planar surfaces
![Generic Plane Benchmark](figures/BM_INTERSECT_PLANES.png)

### mask benchmarks

The following benchmarks test the application of masks on already on-surface intersections.

![Disc Mask Benchmark](figures/BM_DISC2_MASK.png)
![Rectangle Mask Benchmark](figures/BM_RECTANGLE2_MASK.png)
![Ring Mask Benchmark](figures/BM_RING2_MASK.png)
![Trapezoid Mask Benchmark](figures/BM_TRAPEZOID2_MASK.png)


