# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:27:18 2020

"""

import numpy as np
from math import pi,sin,cos
import matplotlib.pyplot as plt
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D


def heat_1d(lx, td, alpha, sig):
    """ This function accepts a linear signal. Then it evolves that signal 
    according to the heat equation. lx: Disccretized length of the domain. tx:
    discretized time."""
    uf = np.tile(sig,(len(td),1))
    
    dt = td[1] - td[0]
    dx = lx[1] - lx[0]    
    r = alpha*dt/dx**2
    r2 = 1-2*r
    
    for t_id,t in enumerate(td[1:]):
        uf[t_id+1,1:-1] = r*uf[t_id,:-2] + r2*uf[t_id,1:-1] + r*uf[t_id,2:]
    
    return uf
    
def heat_2d(lx, td, alpha, sig):
    """ This function accepts a grid signal. Then it evolves that signal 
    according to the heat equation. lx: Disccretized unit of length. tx:
    discretized time."""
    uf = np.repeat(sig[:, :, np.newaxis], len(td), axis=2)
    
    dt = td[1] - td[0]
    dx = lx[1] - lx[0]    
    r = alpha*dt/dx**2
    r2 = 1-4*r
    
    for t_id,_ in enumerate(td[1:]):
        uf[1:-1,1:-1, t_id+1] = r*uf[1:-1,:-2,t_id] + r*uf[1:-1,2:,t_id] + r*uf[2:,1:-1,t_id] + r*uf[:-2,1:-1,t_id] + r2*uf[1:-1,1:-1,t_id]
    
    return uf
    
def heat_2d_save(lx, td, alpha, sig, num):
    """ This function accepts a grid signal. Then it evolves that signal 
    according to the heat equation. lx: Disccretized unit of length. tx:
    discretized time. This function does the same operation as heat_2d but
    saves only a few frames."""
    k = int((len(td)-2)*1.0/(num-1))
    
    uf = np.repeat(sig[:, :, np.newaxis], num, axis=2)
    tmp_im = np.zeros((sig.shape[0],sig.shape[1]))    
    tmp_im[:,:] = sig
    
    dt = td[1] - td[0]
    dx = lx[1] - lx[0]    
    r = alpha*dt/dx**2
    r2 = 1-4*r
    
    for t_id,_ in enumerate(td[1:]):        
        tmp_im[1:-1,1:-1] = r*tmp_im[1:-1,:-2] + r*tmp_im[1:-1,2:] + r*tmp_im[2:,1:-1] + r*tmp_im[:-2,1:-1] + r2*tmp_im[1:-1,1:-1]
        
        rem = t_id%k
        ind = int(t_id/k)
        if rem == 0:
            uf[:,:, ind] = tmp_im
    
    return uf
    
def get_coord(t):
    """ This function gives us the co-ordinates corresponding to the value of 
    the parameter t. t may also be an array."""
    coord = np.zeros(3)
    coord[0] = t*np.cos(6*t)
    coord[1] = t*np.sin(6*t)
    coord[2] = t
    return coord
    
def get_tangent(t):
    """ This function takes the value of t and returns the tangent for the
    curve (tcos6t,tsin6t,t)."""
    tgnt = np.zeros(3)
    tgnt[0] = -6*t*sin(6*t) + cos(6*t)
    tgnt[1] = 6*t*cos(6*t) + sin(6*t)
    tgnt[2] = 1
    
    return tgnt
    
def get_normal(t):
    """ This function takes the value of t and returns the tangent for the
    curve (tcos6t,tsin6t,t)."""
    nrml_tmp = np.zeros(3)
    nrml = np.zeros(3)
    nrml_tmp[0] = -36*t*cos(6*t)-12*sin(6*t)
    nrml_tmp[1] = -36*t*sin(6*t)+12*sin(6*t)
    nrml_tmp[2] = 0
    nrml[:] = nrml_tmp/np.linalg.norm(nrml_tmp)
    
    return nrml
    
def main():
    alpha = 0.2
    L = 10.0
    nL = 100        
    lx = np.linspace(0, L, nL)
    T = 1.0
    nT = 100
    td = np.linspace(0, T, nT)
    sig = np.sin(2*pi*lx)/(lx+0.1) 
    sig[0] = sig[-1] = 0
    # Evolve the signal on the line
    uf = heat_1d(lx, td, alpha, sig)
    # Show the heat evolution wrt time
    plt.figure()
    plt.imshow(uf)
    plt.colorbar()
    plt.savefig('./images/heat1D.png')
    
    # Evolve the signal on the grid
    alpha_2d = 0.2
    T_2d = 0.25
    nT_2d = 100
    sp0 = 0    
    sp1 = int(nT_2d/4.0)
    sp2 = int(nT_2d/2.0)
    sp3 = int(nT_2d)-1
    td_2d = np.linspace(0, T_2d, nT_2d)
    ly = np.linspace(0, L, nL)
    xm, ym = np.meshgrid(lx, ly)
    sig_2d = np.multiply(np.sin(2*pi*xm), np.cos(2*pi*ym))
    sig_2d[:,0] = sig_2d[:,-1] = sig_2d[0,:] = sig_2d[-1,:] = 0
    uf_2d = heat_2d(lx, td_2d, alpha_2d, sig_2d)
    plt.figure()
    plt.imshow(uf_2d[:,:,sp0],vmin=-1.0, vmax=1.0)
    plt.colorbar()
    plt.savefig('./images/heat2D_0.png')
    plt.figure()
    plt.imshow(uf_2d[:,:,sp1],vmin=-1.0, vmax=1.0)
    plt.colorbar()
    plt.savefig('./images/heat2D_1.png')
    plt.figure()
    plt.imshow(uf_2d[:,:,sp2],vmin=-1.0, vmax=1.0)
    plt.colorbar()
    plt.savefig('./images/heat2D_2.png')
    plt.figure()
    plt.imshow(uf_2d[:,:,sp3],vmin=-1.0, vmax=1.0)
    plt.colorbar()
    plt.savefig('./images/heat2D_3.png')
    
    # Evolve the image using heat equation
    alpha_2d_im = 0.2
    T_2d_im = 50.0
    nT_2d_im = 10000
    k = 30    
    sp0_im = 0
    sp2_im = int(k/5.0)
    sp3_im = int(k)-1
    td_2d_im = np.linspace(0, T_2d_im, nT_2d_im)
    im_bw = np.array(Image.open("./data/barbara.png"))/256.0
    im_bw[:,0] = im_bw[:,-1] = im_bw[0,:] = im_bw[-1,:] = 0
    im_bw_out = heat_2d_save(lx, td_2d_im, alpha_2d_im, im_bw, k)
    
    for i in range(im_bw_out.shape[2]):
        new_im = Image.fromarray(np.around(255*im_bw_out[:,:,i])).convert('RGB')
        new_im.save("./images/barbara/barbara_heat_"+str(i).zfill(2)+".png")
    
    # Create movie showing the heat change wrt time
    
    # Plot a Frenet frame for a 3d curve    
    z1 = np.linspace(1,6.0,200)    
    x1 = z1*np.cos(6*z1)
    y1 = z1*np.sin(6*z1)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x1, y1, z1)    
    n_arrow = 7
    t_plot = z1[range(0,len(z1),int(len(z1)/n_arrow))]
    
    xt1 = np.zeros(3*len(t_plot))
    yt1 = np.zeros(3*len(t_plot))
    zt1 = np.zeros(3*len(t_plot))
    xa1 = np.zeros(3*len(t_plot))
    ya1 = np.zeros(3*len(t_plot))
    za1 = np.zeros(3*len(t_plot))
    
    # Plot arrows temp code
    n_vs = 3
    
    # Save the vectors for the Frenet frame
    for t_id,ti in enumerate(t_plot):
        xt1[t_id*n_vs:(t_id+1)*n_vs], yt1[t_id*n_vs:(t_id+1)*n_vs], zt1[t_id*n_vs:(t_id+1)*n_vs] = tuple(get_coord(ti))
        xa1[t_id*n_vs], ya1[t_id*n_vs], za1[t_id*n_vs] = tuple(get_tangent(ti))
        xa1[t_id*n_vs+1], ya1[t_id*n_vs+1], za1[t_id*n_vs+1] = tuple(get_normal(ti))
        xa1[t_id*n_vs+2], ya1[t_id*n_vs+2], za1[t_id*n_vs+2] = tuple(np.cross(get_tangent(ti),get_normal(ti)))        
            
    ax.quiver(xt1, yt1, zt1, xa1, ya1, za1, length=1,color='r') # normalize=True,
    plt.savefig("./images/frenet_3d.png")
    
    plt.show()    
    return
    
if __name__== "__main__":
    main()