import awkward
import coffea
#import coffea.nanoevents
#import coffea.nanoevents.methods
import numba
import numpy
# import awkward.numba

from coffea.nanoevents.methods import candidate
awkward.behavior.update(candidate.behavior)

"""
CaloPoint and CaloPoint3D implementation ported from:
https://github.com/cms-sw/cmssw/blob/master/DataFormats/JetReco/src/Jet.cc
"""


spec_CaloPoint3D = [
    ("vertex", numba.float64[:]),
    ("direction", numba.float64[:]),
    ("point", numba.float64[:]),
    ("depth", numba.float64),
    ("R_BARREL", numba.float64),
    ("R_BARREL2", numba.float64),
    ("Z_ENDCAP", numba.float64),
    ("R_FORWARD", numba.float64),
    ("R_FORWARD2", numba.float64),
    ("Z_FORWARD", numba.float64),
    ("Z_BIG", numba.float64),
]

@numba.experimental.jitclass(spec_CaloPoint3D)
class CaloPoint3D :
    
    def __init__(
        self,
        vertex,
        direction,
    ) :
        
        self.depth = 0.1
        self.R_BARREL = (1. - self.depth) * 143. + self.depth * 407.
        self.R_BARREL2 = self.R_BARREL * self.R_BARREL
        self.Z_ENDCAP = (1. - self.depth) * 320. + self.depth * 568.                         # 1/2(EEz+HEz)
        self.R_FORWARD = self.Z_ENDCAP / numpy.sqrt(numpy.cosh(3.) * numpy.cosh(3.) - 1.)  # eta=3
        self.R_FORWARD2 = self.R_FORWARD * self.R_FORWARD
        self.Z_FORWARD = 1100. + self.depth * 165.
        self.Z_BIG = 1.e5
        
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
        
        
        self.point = numpy.zeros(3)
        self.point[0] = dirUnit[0] * r - dirUnit[1] * tIP
        self.point[1] = dirUnit[1] * r + dirUnit[0] * tIP
        self.point[2] = z


# @numba.njit
# def getMagP3(v) :
    
#     mag = numpy.sqrt(numpy.sum(v**2))
    
#     return mag


@numba.njit
def physicsP4(
    inParticle, # (px, py, pz, E)
    # oldVertex, # (x, y, z)
    newVertex, # (x, y, z)
) :
    """
    Argument format:
        newVertex       -> (x, y, z)
        inParticle      -> (px, py, pz, E)
        oldVertex       -> (x, y, z)
    Description:
        Will change the vertex of "inParticle" from "oldVertex" to "newVertex" and return the new p4(px, py, pz, E).
    """
    

    # inParticle = numpy.array([_val for _val in inParticle])#, dtype = float)
    # oldVertex = numpy.array([_val for _val in oldVertex])#, dtype = float)
    oldVertex = numpy.array([0.0, 0.0, 0.0])
    # newVertex = numpy.array([_val for _val in newVertex])#, dtype = float)
    
    # inParticle_p = inParticle[0: 3].copy()
    
    # Jet position in Calo
    caloPoint = CaloPoint3D(
        #vertex = oldVertex,
        #direction = inParticle_p,
        oldVertex,
        inParticle[0: 3],
    )
    
    physicsDir = caloPoint.point - newVertex

    # p = getMagP3(inParticle_p)
    # physicsDir_unit = physicsDir / getMagP3(physicsDir)
    # p3 = p * physicsDir_unit
    
    # resultP4 = numpy.array([
    #     p3[0],
    #     p3[1],
    #     p3[2],
    #     inParticle[3],
    # ])
    
    # return resultP4

    #temporary return only eta phi at the moment
    p3 = physicsDir
    pt = numpy.sqrt(p3[0]*p3[0]+p3[1]*p3[1])

    return numpy.array([
        numpy.arcsinh(p3[2]/pt),
        numpy.arcsin(p3[1]/pt)
    ])


vphysicsP4 = numpy.vectorize(
    pyfunc = physicsP4,
    signature = "(4),(3)->(2)",
)


@numba.njit
def np_delta_r2(eta_1, phi_1, eta_2, phi_2):
    return numpy.sqrt(numpy.square(eta_1-eta_2)+numpy.square(phi_1-phi_2))


def coffea_nearest_metric_deltaR_shiftVertex(
    v1s,
    v2s,
) :
    
    offsets_v2 = v2s.layout.offsets
    offsets_v2_sub = v2s.layout.content.offsets
    v2_x = numpy.array(v2s.x.layout.content.content)
    v2_y = numpy.array(v2s.y.layout.content.content)
    v2_z = numpy.array(v2s.z.layout.content.content)
    v2_t = numpy.array(v2s.t.layout.content.content)
    v2_momentum = numpy.column_stack((v2_x, v2_y, v2_z, v2_t))

    v1_eta = numpy.array(v1s.eta.layout.content.content)
    v1_phi = numpy.array(v1s.phi.layout.content.content)

    # offsets_v1 = v1s.layout.offsets
    # offsets_v1_sub = v1s.layout.content.offsets
    v1_vx = numpy.array(v1s.vertexX.layout.content.content)
    v1_vy = numpy.array(v1s.vertexY.layout.content.content)
    v1_vz = numpy.array(v1s.vertexZ.layout.content.content)
    v1_vtx = numpy.column_stack((v1_vx, v1_vy, v1_vz))

    v2_shifted = vphysicsP4(
            inParticle = v2_momentum,
            newVertex = v1_vtx,
        )

    result_np = np_delta_r2(v1_eta, v1_phi, v2_shifted[:,0], v2_shifted[:,1])

    results =  awkward.Array(
        awkward.layout.ListOffsetArray64(
            awkward.layout.Index64(offsets_v2),
            awkward.layout.ListOffsetArray64(
                awkward.layout.Index64(offsets_v2_sub),
                awkward.layout.NumpyArray(result_np)
            )
        )
    )

    return results

