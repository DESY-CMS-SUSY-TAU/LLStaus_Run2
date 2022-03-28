import numpy


"""
CaloPoint implementation ported from:
https://github.com/cms-sw/cmssw/blob/master/DataFormats/JetReco/src/Jet.cc
"""

class CaloPoint :
    
    depth = 0.1
    R_BARREL = (1. - depth) * 143. + depth * 407.
    R_BARREL2 = R_BARREL * R_BARREL
    Z_ENDCAP = (1. - depth) * 320. + depth * 568.                         # 1/2(EEz+HEz)
    R_FORWARD = Z_ENDCAP / numpy.sqrt(numpy.cosh(3.) * numpy.cosh(3.) - 1.)  # eta=3
    R_FORWARD2 = R_FORWARD * R_FORWARD
    Z_FORWARD = 1100. + depth * 165.
    Z_BIG = 1.e5


#new implementation to derive CaloPoint for free 3d vertex.
#code provided thanks to Christophe Saout
class CaloPoint3D(CaloPoint) :
    
    def __init__(
        self,
        vertex,
        direction,
    ) :
        # note: no sanity checks here, make sure vertex is inside the detector!
    
        # check if positive or negative (or none) endcap should be tested
        side = None
        
        if (direction[2] < -1e-9) :
            side = -1
        elif (direction[2] > 1e-9) :
            side = +1
        else :
            side = 0
        
        
        dirR = (direction[0]**2.0 + direction[1]**2.0)**0.5
    
        # normalized direction in x-y plane
        dirUnit = (direction[0] / dirR, direction[1] / dirR)
    
        # rotate the vertex into a coordinate system where direction lies along x
    
        # vtxLong is the longitudinal coordinate of the vertex wrt/ direction
        vtxLong = dirUnit[0] * vertex[0] + dirUnit[1] * vertex[1]
    
        # tIP is the (signed) transverse impact parameter
        tIP = dirUnit[0] * vertex[1] - dirUnit[1] * vertex[0]
    
        # r and z coordinate
        r = None
        z = None
    
        if (side) :
            slope = dirR / direction[2]
    
        # check extrapolation to endcap
        r = vtxLong + slope * (side * self.Z_ENDCAP - vertex[2])
        r2 = r**2 + tIP**2
    
        if (r2 < self.R_FORWARD2) :
            # we are in the forward calorimeter, recompute
            r = vtxLong + slope * (side * self.Z_FORWARD - vertex[2])
            z = side * self.Z_FORWARD
        elif (r2 < self.R_BARREL2) :
            # we are in the endcap
            z = side * self.Z_ENDCAP
        else :
            # we are in the barrel, do the intersection below
            side = 0
    
        if (not side) :
            # we are in the barrel
            slope = direction[2] / dirR
            r = numpy.sqrt(self.R_BARREL2 - tIP**2)
            z = vertex[2] + slope * (r - vtxLong)
        
    
        # rotate (r, tIP, z) back into original x-y coordinate system
        self.point = numpy.array([
            dirUnit[0] * r - dirUnit[1] * tIP,
            dirUnit[1] * r + dirUnit[0] * tIP,
            z
        ])
    

def getMagP3(v) :
    
    mag = numpy.sqrt(numpy.sum(v**2))
    
    return mag


def physicsP4(
    newVertex, # (x, y, z)
    inParticle, # (px, py, pz, E)
    oldVertex = (0.0, 0.0, 0.0), # (x, y, z)
) :
    """
    Argument format:
        newVertex       -> (x, y, z)
        inParticle      -> (px, py, pz, E)
        oldVertex       -> (x, y, z)
    Description:
        Will change the vertex of "inParticle" from "oldVertex" to "newVertex" and return the new p4(px, py, pz, E).
    """
    
    newVertex = numpy.array(newVertex, dtype = float)
    inParticle = numpy.array(inParticle, dtype = float)
    oldVertex = numpy.array(oldVertex, dtype = float)
    
    inParticle_p = numpy.array(inParticle[0: 3])
    
    # Jet position in Calo
    caloPoint = CaloPoint3D(
        vertex = oldVertex,
        direction = inParticle_p,
    )
    
    physicsDir = caloPoint.point - newVertex
    p = getMagP3(inParticle_p)
    physicsDir_unit = physicsDir / getMagP3(physicsDir)
    p3 = p * physicsDir_unit
    
    physicsP4 = numpy.array([
        p3[0],
        p3[1],
        p3[2],
        inParticle[3],
    ])
    
    return physicsP4

