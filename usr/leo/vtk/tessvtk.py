import math
import numpy
from mayavi import mlab
from tvtk.api import tvtk


def sph2cart(lon,lat,r):
    R = 6378137.0
    d2r = math.pi/180.0
    r2d = 180.0/math.pi
    x = math.cos(d2r*lat)*math.cos(d2r*lon)*(R + r)
    y = math.cos(d2r*lat)*math.sin(d2r*lon)*(R + r)
    z = math.sin(d2r*lat)*(R + r)
    return x, y, z

def tess2vtk(w,e,s,n,t,b):
    points = [sph2cart(w,s,b),
              sph2cart(e,s,b),
              sph2cart(e,n,b),
              sph2cart(w,n,b),
              sph2cart(w,s,t),
              sph2cart(e,s,t),
              sph2cart(e,n,t),
              sph2cart(w,n,t),
              sph2cart(0.5*(w+e),s,b),
              sph2cart(e,0.5*(s+n),b),
              sph2cart(0.5*(w+e),n,b),
              sph2cart(w,0.5*(s+n),b),
              sph2cart(0.5*(w+e),s,t),
              sph2cart(e,0.5*(s+n),t),
              sph2cart(0.5*(w+e),n,t),
              sph2cart(w,0.5*(s+n),t),
              sph2cart(w,s,0.5*(t+b)),
              sph2cart(e,s,0.5*(t+b)),
              sph2cart(e,n,0.5*(t+b)),
              sph2cart(w,n,0.5*(t+b)),
             ]
    for p in points:
        print p[0], p[1], p[2]
    cells = [20]
    cells.extend(range(20))
    offsets = [0]
    scalars = [1]
    cell_array = tvtk.CellArray()
    cell_array.set_cells(1, numpy.array(cells))
    cell_types = numpy.array([25])
    vtkmesh = tvtk.UnstructuredGrid(points=numpy.array(points))
    vtkmesh.set_cells(cell_types, numpy.array(offsets), cell_array)
    vtkmesh.cell_data.scalars = numpy.array(scalars)
    return vtkmesh

if __name__ == '__main__':

    tess = [30,60,-15,15,0,-3330000]
    mesh = tess2vtk(*tess)
    mlab.pipeline.surface(mesh)
    mlab.show()

