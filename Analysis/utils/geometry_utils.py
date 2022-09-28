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

######################
######################
# non jit version below
######################
######################
'''
import awkward
import coffea
#import coffea.nanoevents
#import coffea.nanoevents.methods
import numba
import numpy

from coffea.nanoevents.methods import candidate
awkward.behavior.update(candidate.behavior)

"""
CaloPoint and CaloPoint3D implementation ported from:
https://github.com/cms-sw/cmssw/blob/master/DataFormats/JetReco/src/Jet.cc
"""

#@numba.jitclass
#class CaloPoint :
#    
#    depth = 0.1
#    R_BARREL = (1. - depth) * 143. + depth * 407.
#    R_BARREL2 = R_BARREL * R_BARREL
#    Z_ENDCAP = (1. - depth) * 320. + depth * 568.                         # 1/2(EEz+HEz)
#    R_FORWARD = Z_ENDCAP / numpy.sqrt(numpy.cosh(3.) * numpy.cosh(3.) - 1.)  # eta=3
#    R_FORWARD2 = R_FORWARD * R_FORWARD
#    Z_FORWARD = 1100. + depth * 165.
#    Z_BIG = 1.e5


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
class CaloPoint3D : #(CaloPoint) :
    
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
        
    
        # rotate (r, tIP, z) back into original x-y coordinate system
        #self.point = numpy.array([
        #    dirUnit[0] * r - dirUnit[1] * tIP,
        #    dirUnit[1] * r + dirUnit[0] * tIP,
        #    z,
        #], dtype = numpy.float)
        
        self.point = numpy.zeros(3)
        self.point[0] = dirUnit[0] * r - dirUnit[1] * tIP
        self.point[1] = dirUnit[1] * r + dirUnit[0] * tIP
        self.point[2] = z


@numba.njit
def getMagP3(v) :
    
    mag = numpy.sqrt(numpy.sum(v**2))
    
    return mag


@numba.njit
def physicsP4(
    inParticle, # (px, py, pz, E)
    oldVertex, # (x, y, z)
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
    
    #arr = numpy.array(numba.typed.List([0.8, 0.6, -1.7]))
    
    #print("inParticle:", inParticle)
    #print("oldVertex:", oldVertex)
    #print("newVertex:", newVertex)
    
    #inParticle = numpy.array(inParticle, dtype = float)
    #oldVertex = numpy.array(oldVertex, dtype = float)
    #newVertex = numpy.array(newVertex, dtype = float)
    
    inParticle = numpy.array([_val for _val in inParticle])#, dtype = float)
    oldVertex = numpy.array([_val for _val in oldVertex])#, dtype = float)
    newVertex = numpy.array([_val for _val in newVertex])#, dtype = float)
    
    ##inParticle_p = numpy.array(inParticle[0: 3])
    inParticle_p = inParticle[0: 3].copy()
    
    # Jet position in Calo
    caloPoint = CaloPoint3D(
        #vertex = oldVertex,
        #direction = inParticle_p,
        oldVertex,
        inParticle_p,
    )
    
    physicsDir = caloPoint.point - newVertex
    p = getMagP3(inParticle_p)
    physicsDir_unit = physicsDir / getMagP3(physicsDir)
    p3 = p * physicsDir_unit
    
    resultP4 = numpy.array([
        p3[0],
        p3[1],
        p3[2],
        inParticle[3],
    ])
    
    return resultP4


vphysicsP4 = numpy.vectorize(
    pyfunc = physicsP4,
    signature = "(4),(3),(3)->(4)",
    #signature = "(4),(3),(nvtx,3)->(nvtx,4)",
    #excluded = ["newVertex", "oldVertex"],
)


def numpy_to_coffeaLorentzVector(x, y, z, t) :
    
    lv = awkward.zip(
        {
            "x": x,
            "y": y,
            "z": z,
            "t": t,
        },
        with_name = "LorentzVector"
        #parameters = {"__record__": "LorentzVector"}
    )
    
    return lv


def get_p4_shiftVertex(
    v1s,
    v2s,
    oldVertex = None,
    newVertex = None,
) :
    """
    Change the vertex of v2s from "oldVertex" to "newVertex" and return the deltaR w.r.t. v1s.
    Will attempt to use (vertexX, vertexY, vertexZ) from v2s if oldVertex is not provided.
    Will attempt to use (vertexX, vertexY, vertexZ) from v1s if newVertex is not provided.
    Designed to work with the coffea "nearest" metric which broadcasts v1 and v2 such that,
    v1s = [
        [
            [v1_0, v1_0, v1_0],
            [v1_1, v1_1, v1_1]
        ], # Event 0
        ...
    ]
    
    v2s = [
        [
            [v2_0, v2_1, v2_2],
            [v2_0, v2_1, v2_2],
        ], # Event 0
        ...
    ]
    
    when Event 0 has 2 entries for v1 and 3 entries for v2.
    """
    
    builder = awkward.ArrayBuilder()
    
    assert(len(v1s) == len(v2s))
    
    for iEvt in range(len(v1s)) :
        
        with builder.list() :
            
            vv1s = v1s[iEvt]
            vv2s = v2s[iEvt]
            
            assert(len(vv1s) == len(vv2s))
            
            for iv in range(len(vv1s)) :
                
                vvv1s = vv1s[iv]
                vvv2s = vv2s[iv]
                
                #print("vvv1s:", vvv1s)
                #print("vvv2s:", vvv2s)
                
                assert(len(vvv1s) == len(vvv2s))
                
                v2_shifted = []
                
                if (len(vvv2s)) :
                    
                    v2_in = awkward.to_list(awkward.zip([vvv2s.x, vvv2s.y, vvv2s.z, vvv2s.t]))
                    oldVtxs = oldVertex
                    newVtxs = newVertex
                    
                    if (oldVtxs is None) :
                        
                        oldVtxs = awkward.to_list(awkward.zip([vvv2s.vertexX, vvv2s.vertexY, vvv2s.vertexZ]))
                    
                    if (newVtxs is None) :
                        
                        newVtxs = awkward.to_list(awkward.zip([vvv1s.vertexX, vvv1s.vertexY, vvv1s.vertexZ]))
                    
                    v2_shifted = vphysicsP4(
                        inParticle = v2_in,
                        oldVertex = oldVertex,
                        newVertex = newVtxs,
                    )
                
                with builder.list() :
                    
                    for v2_tmp in v2_shifted :
                        
                        builder.begin_record("LorentzVector")
                        builder.field("x").real(v2_tmp[0])
                        builder.field("y").real(v2_tmp[1])
                        builder.field("z").real(v2_tmp[2])
                        builder.field("t").real(v2_tmp[3])
                        builder.end_record()
    
    result = builder.snapshot()
    
    return result


def coffea_nearest_metric_deltaR_shiftVertex(
    v1s,
    v2s,
    oldVertex = None,
    newVertex = None,
) :
    """
    Designed to work as a 
    """
    
    v2s_shifted = get_p4_shiftVertex(
        v1s,
        v2s,
        oldVertex,
        newVertex,
    )
    
    result = v1s.delta_r(v2s_shifted)
    #print("delta_r:", result)
    
    return result


#def get_deltaR_shiftVertex(
#    v1s,
#    v2s,
#    oldVertex = None,
#    newVertex = None,
#) :
#    """
#    Change the vertex of v2s from "oldVertex" to "newVertex" and return the deltaR w.r.t. v1s.
#    Will attempt to use (vertexX, vertexY, vertexZ) from v2s if oldVertex is not provided.
#    Will attempt to use (vertexX, vertexY, vertexZ) from v1s if newVertex is not provided.
#    """
#    
#    builder = awkward.ArrayBuilder()
#    
#    assert(len(v1s) == len(v2s))
#    
#    for iEvt in range(len(v1s)) :
#        
#        builder.begin_list()
#        
#        vv1s = v1s[iEvt]
#        vv2s = v2s[iEvt]
#        
#        assert(len(vv1s) == len(vv2s))
#        
#        for iv in range(len(vv1s)) :
#            
#            vvv1s = vv1s[iv]
#            vvv2s = vv2s[iv]
#            
#            assert(len(vvv1s) == len(vvv2s))
#            
#            #print("vvv1s:", vvv1s)
#            #print("vvv2s:", vvv2s)
#            #print(len(vvv1s), len(vvv2s))
#            
#            v2_in = awkward.to_list(awkward.zip([vvv2s.x, vvv2s.y, vvv2s.z, vvv2s.t]))
#            oldVtxs = oldVertex
#            newVtxs = newVertex
#            
#            if (oldVtxs is None) :
#                
#                oldVtxs = awkward.to_list(awkward.zip([vvv2s.vertexX, vvv2s.vertexY, vvv2s.vertexZ]))
#            
#            if (newVtxs is None) :
#                
#                newVtxs = awkward.to_list(awkward.zip([vvv1s.vertexX, vvv1s.vertexY, vvv1s.vertexZ]))
#            
#            print("v2_in:", v2_in)
#            #print("newVtxs:", newVtxs)
#            
#            v2_shifted = vphysicsP4(
#                inParticle = v2_in,
#                oldVertex = oldVertex,
#                newVertex = newVtxs,
#            )
#            
#            print("v2_shifted:", v2_shifted)
#            
#            v2_shiftedP4s = awkward.zip(
#                {
#                    "x": v2_shifted[:, 0],
#                    "y": v2_shifted[:, 1],
#                    "z": v2_shifted[:, 2],
#                    "t": v2_shifted[:, 3],
#                },
#                with_name = "LorentzVector"
#            )
#            
#            #print("v2_shiftedP4s:", v2_shiftedP4s)
#            
#            arr_dR = vvv1s.delta_r(v2_shiftedP4s)
#            #print("arr_dR:", arr_dR)
#            
#            builder.begin_list()
#            
#            for dR in arr_dR :
#                
#                builder.real(dR)
#            
#            builder.end_list()
#        
#        builder.end_list()
#    
#    result = builder.snapshot()
#    print("result:", result)
#    
#    return result
'''