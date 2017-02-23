# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 08:58:21 2017

@author: Callum
"""

"""
Example script shoing conversion of a dxf file to a poly object, saving the 
object to a .poly file, and displaying the PSLG
"""

from poly import poly
from polyGui import polyGui

PSLG = poly(filePath = "Demo.dxf",elementSize = 2e-5,writeFile=True)
pGui = polyGui(PSLG,title='Converted data from .dxf file')

PSLGFromPoly = poly(filePath = "Demo.dxf.poly",elementSize = 2e-5,writeFile=False)
pGuiFromPoly = polyGui(PSLGFromPoly,title='Read data from .poly file')

PSLGFromPoly2 = poly(filePath = "Demo2.poly",elementSize = 2e-5,writeFile=False)
pGuiFromPoly2 = polyGui(PSLGFromPoly2,title='Read data from .poly file')





