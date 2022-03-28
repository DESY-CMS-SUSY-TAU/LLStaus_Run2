import geometry_utils
import ROOT


newVertex = (20.0, 20.0, 0.0)
detP4_PtEtaPhiM = (50.0, 0.7, 1.9, 10.4)
oldVertex = (0.0, 0.0, 0.0)

rt_detP4 = ROOT.Math.PtEtaPhiMVector(*detP4_PtEtaPhiM)

detP4 = (rt_detP4.Px(), rt_detP4.Py(), rt_detP4.Pz(), rt_detP4.E())

physP4 = geometry_utils.physicsP4(
    newVertex = newVertex,
    inParticle = detP4,
    oldVertex = oldVertex,
)

rt_phP4 = ROOT.Math.PxPyPzEVector(*physP4)


print("Old: pt %0.4f, eta %+0.4f, phi %+0.4f, mass %0.4f, energy %0.4f" %(rt_detP4.Pt(), rt_detP4.Eta(), rt_detP4.Phi(), rt_detP4.M(), rt_detP4.E()))
print("New: pt %0.4f, eta %+0.4f, phi %+0.4f, mass %0.4f, energy %0.4f" %(rt_phP4.Pt(), rt_phP4.Eta(), rt_phP4.Phi(), rt_phP4.M(), rt_phP4.E()))
