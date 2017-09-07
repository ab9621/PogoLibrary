import numpy as np

def ThreePointToCenterAndAngles(p1,p2,p3):
    def _getDiffs():
        dy1 = float(p2[1]-p1[1])
        dx1 = float(p2[0]-p1[0])
        dy2 = float(p3[1]-p2[1])
        dx2 = float(p3[0]-p2[0])
        return dy1,dx1,dy2,dx2
    dy1,dx1,dy2,dx2 = _getDiffs()
    if dx1 == 0:
        p2,p3 = p3,p2
        dy1,dx1,dy2,dx2 = _getDiffs()
    elif dx2 == 0:
        p1,p2 = p2,p1
        dy1,dx1,dy2,dx2 = _getDiffs()
        
    grad1 = dy1/dx1
    grad2 = dy2/dx2
    denom = 2*(grad2-grad1)
    if denom == 0 or np.any([dx1,dx2]==0):
        errMsg = 'Points lie on parallel lines. No circle could be found'
        raise ValueError(errMsg)
    numerator = grad1*grad2*(p1[1]-p3[1])+grad2*(p1[0]+p2[0])-grad1*(p2[0]+p3[0])
    x = numerator / denom
    y = -1/grad1*(x-(p1[0]+p2[0])/2)+(p1[1]+p2[1])/2
    
    
    th1,th2,th3 = getAnglesRelativeToCenter(p1,[x,y]),getAnglesRelativeToCenter(p2,[x,y]),getAnglesRelativeToCenter(p3,[x,y])
    startTheta = min([th1,th2,th3])
    endTheta = max([th1,th2,th3])
    
    radius = np.sqrt((p1[0]-x)*(p1[0]-x)+(p1[1]-y)*(p1[1]-y))
    return x,y,startTheta,endTheta,radius

def SCEToCentreAndAngles(p1,p2,p3):
    startTheta,endTheta = getAnglesRelativeToCentre(p1,p2),getAnglesRelativeToCentre(p3,p2)
    radius = np.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]))
    return p2[0],p2[1],startTheta,endTheta,radius

def SCAToCentreAndAngles(p1,p2,p3):
    startTheta = getAnglesRelativeToCentre(p1,p2)
    radius = np.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]))
    return p2[0],p2[1],startTheta,p3,radius
    
def SCLToCentreAndAngles(p1,p2,p3):
    startTheta = getAnglesRelativeToCentre(p1,p2)
    radius = np.sqrt((p1[0]-p2[0])**2+(p1[1]+p2[1])**2)
    dTh = p3/radius
    radius = np.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]))
    return p2[0],p2[1],startTheta,startTheta+dTh,radius

def getAnglesRelativeToCenter(point,centre):
    dy = point[1]-centre[1]    
    dx = point[0]-centre[0]
    th = np.arctan2(dy,dx)
    return th