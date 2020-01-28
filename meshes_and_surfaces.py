# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:54:52 2020

"""

from skimage import measure
from mayavi import mlab
import visvis as vv
from skimage.draw import ellipsoid
from math import pi
import numpy as np
from read_off import read_off
import sys
import math

def triangle_compute_area(v1,v2,v3):
    """This function computes the area of a triangle given the co-ordinates of 
    its three vertices."""
    
    return np.linalg.norm(np.cross(v2-v1,v3-v1))/2.0
    
def triangle_get_normal(v1, v2, v3):
    """ This function computes the normal for a triangle given by v1, v2, v3.
    Then it returns the center of the triangle and the normal of the triangle."""
    
    
    return np.mean(np.array([v1, v2, v3]), axis=0), np.cross(v2-v1,v3-v1)

def mesh_get_normals(vs, fs):
    if not isinstance(fs,np.ndarray) or not isinstance(vs,np.ndarray):
        raise Exception('The vertices and faces should be N,3 numpy arrays.')
        
    n_faces = fs.shape[0]
    centers = np.zeros((n_faces,3))
    normals = np.zeros((n_faces,3))
    
    for i,face in enumerate(fs):
        centers[i,:], normals[i,:] = triangle_get_normal(vs[face[0]],vs[face[1]],vs[face[2]])
    
    return centers, normals

def mesh_comp_area(verts,faces):
    area = 0

    # Check the format of vertices and faces
    if not (verts.shape[1] == 3 or faces.shape[1] == 3):
        raise IOError('Size of vertices or faces array not correct')
    
    for face in faces:
        area = area + triangle_compute_area(verts[face[0],:], verts[face[1],:], verts[face[2],:])
    return area    
    
def rot_z(v,theta):
    """This function rotates a numpy array or 3 co-ordinates around the z-axis."""
    
    # Check whether v is numpy array and is 3d
    if (len(v) != 3) or (not isinstance(v,np.ndarray)):
        raise Exception('The vector is not numpy array')
    ct = math.cos(theta)
    st = math.sin(theta)
    rot_mat = np.array([[ct,-st,0],[st,ct,0],[0,0,1]])
    
    return np.dot(rot_mat,v)
    
def check_inside_cylinder(h,r,v):
    """ This function checks if the point v is inside a cylender of height h and
    radius r centered at the origin. Returns 1 if inside and -1 if outside. The
    cylinder has its height along the z-axis."""
    is_inside = -1
    if v[2] < h/2 and v[2] > -h/2 and v[0]**2+v[1]**2<r**2:
        is_inside = 1
    return is_inside
    
def cylinder_iso_fn(h,r,v):
    surf1 = h/2.0-v[2] 
    surf2 = v[2]+h/2.0
    surf3 = r**2-v[0]**2-v[1]**2
    sign = min(np.sign(surf1),np.sign(surf2),np.sign(surf3))
    val = sign*abs(surf1*surf2*surf3)
    return val
    
def torus_iso_fn(R,r,v):
    val = r**2 - v[2]**2 - (math.sqrt(v[0]**2+v[1]**2)-R)**2
    return val   
    
def count_edges(faces):
    """ This function counts the number of edges in a mesh when the datastructure 
    for faces is passed an argument."""
    if not faces.shape[1] == 3:
        raise Exception('We expect a (N,3) size array')
    edge_set = set()
    for face in faces:
        min_v1 = min(face[0],face[1])
        max_v1 = max(face[0],face[1])
        min_v2 = min(face[1],face[2])
        max_v2 = max(face[1],face[2])
        min_v3 = min(face[0],face[2])
        max_v3 = max(face[0],face[2])
        edge_set.add((min_v1,max_v1))
        edge_set.add((min_v2,max_v2))
        edge_set.add((min_v3,max_v3))
    return len(edge_set)

def main():
    # Generate Sphere
    rad = 6
    vertsS, facesS, _, _ = measure.marching_cubes_lewiner(ellipsoid(rad, rad, rad, levelset=True), 0.0)
    mlab.triangular_mesh([vert[0] for vert in vertsS],
                         [vert[1] for vert in vertsS],
                        [vert[2] for vert in vertsS],
                        facesS, colormap='Blues')    
    mlab.savefig('./images/hw1Sphere.png')
    mlab.show()    
    
    # Compute normal for each triangle
    centersS, normalsS = mesh_get_normals(vertsS, facesS)
    # Display normals
    mlab.triangular_mesh([vert[0] for vert in vertsS],
                         [vert[1] for vert in vertsS],
                        [vert[2] for vert in vertsS],
                        facesS, colormap='Blues')
    mlab.quiver3d(centersS[:,0],centersS[:,1],centersS[:,2],normalsS[:,0],normalsS[:,1],normalsS[:,2])
    mlab.savefig('./images/hw1SphereNormals.png')
    mlab.show()
               
    # Number of triangles, faces, unique edges, Euler number
    n_faces_S = facesS.shape[0]
    n_vertices_S = vertsS.shape[0]
    n_edges_S = count_edges(facesS)
    print('The sphere has {0}, {1}, {2} number of triangles, vertices and edges.'.format(n_faces_S, n_vertices_S, n_edges_S))
    print("The sphere's Euler characteristic is {0}.".format(n_vertices_S-n_edges_S+n_faces_S))
    
    
    # Compute surface area
    print('Library computed surface area of Sphere is {0}.'.format(measure.mesh_surface_area(vertsS, facesS)))
    print('Theoretically computed surface area of Sphere is {0}.'.format(4*pi*rad**2))
    print('Mesh computed computed surface area of Sphere is {0}.'.format(mesh_comp_area(vertsS,facesS)))
        
    
    
    # Generate Cylinder
    cy_ht = 10.0
    cy_rad = 1.0
    xr = np.linspace(-1.5,1.5,20)
    yr = np.linspace(-1.5,1.5,20)
    zr = np.linspace(-6,6,35)
    cylinder_vol = np.zeros((len(xr),len(yr),len(zr)))

    for i,xi in enumerate(xr):
        for j,yj in enumerate(yr):
            for k,zk in enumerate(zr):
                cylinder_vol[i,j,k] = cylinder_iso_fn(cy_ht,cy_rad,np.array([xi,yj,zk]))

    vertsC, facesC, _, _ = measure.marching_cubes_lewiner(cylinder_vol, 0.0,spacing=(xr[1]-xr[0],yr[1]-yr[0],zr[1]-zr[0]))
    mlab.triangular_mesh([vert[0] for vert in vertsC],
                         [vert[1] for vert in vertsC],
                        [vert[2] for vert in vertsC],
                        facesC, colormap='Blues')
    mlab.savefig('./images/hw1Cylinder.png')
    mlab.show()        
    
    # Compute normal for each triangle
    centersC, normalsC = mesh_get_normals(vertsC, facesC)
    # Display normals
    mlab.triangular_mesh([vert[0] for vert in vertsC],
                         [vert[1] for vert in vertsC],
                        [vert[2] for vert in vertsC],
                        facesC, colormap='Blues')
    mlab.quiver3d(centersC[:,0],centersC[:,1],centersC[:,2],-normalsC[:,0],-normalsC[:,1],-normalsC[:,2])
    mlab.savefig('./images/hw1CylinderNormals.png')
    mlab.show()
    
    # Count Num triangles, edges, vertices
    n_faces_C = facesC.shape[0]
    n_vertices_C = vertsC.shape[0]
    n_edges_C = count_edges(facesC)
    print('The cylinder has {0}, {1}, {2} number of triangles, vertices and edges.'.format(n_faces_C, n_vertices_C, n_edges_C))
    print("The cylinder's Euler characteristic is {0}.".format(n_vertices_C-n_edges_C+n_faces_C))
    
    # Compute surface area
    print('Library computed surface area of Cylinder is {0}.'.format(measure.mesh_surface_area(vertsC, facesC)))
    print('Theoretically computed surface area of Cylinder is {0}.'.format(2*pi*cy_rad*cy_ht+2*pi*cy_rad**2))
    print('Mesh computed computed surface area of Cylinder is {0}.'.format(mesh_comp_area(vertsC,facesC)))                        
    
    # Import and display surface mesh model from surface repository
    sharp_f = open('./data/104_sharp_sphere.off')
    sharp_vs, sharp_fs = read_off(sharp_f)
    sharp_f.close()
    mlab.triangular_mesh([vert[0] for vert in sharp_vs],
                         [vert[1] for vert in sharp_vs],
                        [vert[2] for vert in sharp_vs],
                        sharp_fs, colormap='Blues')
    mlab.axes()
    mlab.savefig('./images/hw1Dataset.png')
    mlab.show() 

    # Number of triangles, faces, unique edges, Euler number
    n_faces_SS = len(sharp_fs)
    n_vertices_SS = len(sharp_vs)
    n_edges_SS = count_edges(np.array(sharp_fs))
    print('The dataset model has {0}, {1}, {2} number of triangles, vertices and edges.'.format(n_faces_SS, n_vertices_SS, n_edges_SS))
    print("The dataset model's Euler characteristic is {0}.".format(n_vertices_SS-n_edges_SS+n_faces_SS))
                        
    # Rotate by 45 degrees about z-axis
    theta = math.pi/4
    sharp_vs_rot = [list(rot_z(np.array(v),theta)) for v in sharp_vs]             
    mlab.triangular_mesh([vert[0] for vert in sharp_vs_rot],
                         [vert[1] for vert in sharp_vs_rot],
                        [vert[2] for vert in sharp_vs_rot],
                        sharp_fs, colormap='Blues')
    mlab.axes()
    mlab.savefig('./images/hw1DatasetRot.png')
    mlab.show()                        
    

    # Apply 2x scaling in x direction
    mlab.triangular_mesh([2*vert[0] for vert in sharp_vs],
                         [vert[1] for vert in sharp_vs],
                        [vert[2] for vert in sharp_vs],
                        sharp_fs, colormap='Blues')
    mlab.axes()
    mlab.savefig('./images/hw1DatasetSca.png')
    mlab.show()    
    

    # Generate torus mesh
    tor_R = 4.0
    tor_r = 1.5
    xr = np.linspace(-6,6,70)
    yr = np.linspace(-6,6,70)
    zr = np.linspace(-2,2,40)
    tor_vol = np.zeros((len(xr),len(yr),len(zr)))

    for i,xi in enumerate(xr):
        for j,yj in enumerate(yr):
            for k,zk in enumerate(zr):
                tor_vol[i,j,k] = torus_iso_fn(tor_R,tor_r,np.array([xi,yj,zk]))

    vertsT, facesT, _, _ = measure.marching_cubes_lewiner(tor_vol, 0.0)
    mlab.triangular_mesh([vert[0] for vert in vertsT],
                         [vert[1] for vert in vertsT],
                        [vert[2] for vert in vertsT],
                        facesT, colormap='Blues')
    mlab.savefig('./images/hw1Torus.png')
    mlab.show()
    
    return

if __name__== "__main__":
    main()