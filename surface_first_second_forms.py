# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 23:32:53 2020

@author: jinkx
"""

import numpy as np
from dfsio import readdfs
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def main():
    # Import the brain map
    ss = readdfs("./harmonic_mapping_code/surf.cortex.dfs")

    # Generate the u,v map
    uv_map = np.loadtxt("./data/uv_brain_map.txt")

    # Interpolate the 3d coordinates on the 2d grid
    n_pts = 100
    xmin = np.min(uv_map[:,0])
    xmax = np.max(uv_map[:,0])
    ymin = np.min(uv_map[:,1])
    ymax = np.max(uv_map[:,1])
    xg,yg = np.meshgrid(np.linspace(xmin,xmax,n_pts),np.linspace(ymin,ymax,n_pts))
    
    # Interpolating the x, y,z coordiantes
    grid_x = griddata(uv_map, ss.vertices[:,0], (xg, yg), method='linear')
    grid_y = griddata(uv_map, ss.vertices[:,1], (xg, yg), method='linear')
    grid_z = griddata(uv_map, ss.vertices[:,2], (xg, yg), method='linear')

    # Compute the differences
    diffx_u = np.diff(grid_x, axis=0)
    diffx_v = np.diff(grid_x, axis=1)
    diffy_u = np.diff(grid_y, axis=0)
    diffy_v = np.diff(grid_y, axis=1)
    diffz_u = np.diff(grid_z, axis=0)
    diffz_v = np.diff(grid_z, axis=1)
    
    # Compute difference in the u and v direction for the 3d coordinates fn
    alpha_u = np.stack((diffx_u[:,:-1], diffy_u[:,:-1], diffz_u[:,:-1]), axis=2)
    alpha_v = np.stack((diffx_v[:-1,:], diffy_v[:-1,:], diffz_v[:-1,:]), axis=2)
    
    # Compute E,F,G for the nodes in the grid
    E = np.squeeze(np.sum(np.multiply(alpha_u, alpha_u), axis=2))
    F = np.squeeze(np.sum(np.multiply(alpha_u, alpha_v), axis=2))
    G = np.squeeze(np.sum(np.multiply(alpha_v, alpha_v), axis=2))

    # Visualize E,F and G for the square
    plt.figure()
    plt.imshow(E)
    plt.colorbar()
    plt.title("E over the square")
    plt.savefig("./images/E_brain.png")
    
    plt.figure()
    plt.imshow(F)
    plt.colorbar()
    plt.title("F over the square")
    plt.savefig("./images/F_brain.png")
    
    plt.figure()
    plt.imshow(G)
    plt.colorbar()
    plt.title("G over the square")
    plt.savefig("./images/G_brain.png")
    
    plt.show()
    
    return

if __name__== "__main__":
    main()