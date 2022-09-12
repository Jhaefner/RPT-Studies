#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 19:17:26 2022

@author: jonah
"""


#This script executes the secand method using two sets of (x,y) coordinates.

from cmath import sqrt

def secant_method(p1,p2,target):
    """
    Inputs:
    p1 = (x,y) pair, first calculated points
    p2 = (x,y) pair, second calcualted points
    target = desired y-value
    
    Outputs:
    x3 = desired next approximation
    """
    x1=p1[0]; x2=p2[0];
    y1 = p1[1]-target; y2=p2[1]-target;
    x3 = x2-y2*((x2-x1)/(y2-y1));
    return x3