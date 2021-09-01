# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:30:26 2019

@author: Radhanath
"""
from pymol import cmd, stored
import anglebetweenhelices
def findAngle():
    prot = ["4EIY","4MBS"]
    resid = [[283,291,293,305],[337,345,349,356]]
    for x in range(0,len(prot)):
        cmd.fetch("%s" % (prot[x]))
        cmd.select("hel1", "(chain A) and (resi %s-$s)" % (resid[x][0],resid[x][1]))
        cmd.select("hel2", "(chain A) and (resi %s-$s)" % (resid[x][2],resid[x][3]))
        print ("calculating angle for protein %s" % (prot[x]))
        angle_between_helices ("hel1","hel2","helix_orientation",1,0)
        cmd.delete("all")
cmd.extend("findangles", findAngle)