# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 19:10:50 2020
This is to test the differential geometry functions.
"""

import meshes_and_surfaces as ms
import numpy as np
import math

def test_triangle_compute_area():
    v1 = np.array([0,0,0])
    v2 = np.array([0,1,0])
    v3 = np.array([1,0,0])
    if not (ms.triangle_compute_area(v1,v2,v3) == 0.5):
        raise Exception('Computed area should be half.')
    return
    
def test_rot_z():
    v = np.array([1,0,0])
    theta = math.pi/2.0
    vy = np.array([0,1,0])
    if not (abs(ms.rot_z(v,theta) - vy) < 1e-6).all():
        raise Exception('Rotated vector not as expected')
    return    
    
def test_check_inside_cylinder():
    v1 = np.array([0.2,0.2,0])
    v2 = np.array([0.2,0.2,20])
    if ms.check_inside_cylinder(10,1,v1) != 1:
        raise Exception('This point is inside the cylinder')
    if ms.check_inside_cylinder(10,1,v2) != -1:
        raise Exception('This point is outside the cylinder')        
    return    
    
def test_cylinder_iso_fn():
    v1 = np.array([0.2,0.2,0])
    v2 = np.array([0.2,0.2,20])
    if ms.check_inside_cylinder(10,1,v1) < 0:
        raise Exception('This point is inside the cylinder')
    if ms.check_inside_cylinder(10,1,v2) > 0:
        raise Exception('This point is outside the cylinder')        
    return
    
def test_count_edges():
    f1 = np.array([[1,2,3]])
    if ms.count_edges(f1) != 3:
        raise Exception('A simple triangle should have 3 edges')
    return    
    
def main():
#    test_triangle_compute_area()
#    test_rot_z()
#    test_check_inside_cylinder()
#    test_cylinder_iso_fn()
    test_count_edges()
    return
    
if __name__== "__main__":
    main()    