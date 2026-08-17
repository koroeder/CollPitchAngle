import numpy as np

forces = ["f0","f1","f2","f3","f4","f5","f6"]

kT = 0.596

for force in forces:
    data = np.genfromtxt("rho_{}.dat".format(force))
    energy = np.genfromtxt("e_{}.dat".format(force))
    minE = np.min(energy)
    energy = energy - minE
    weights = np.exp(-energy/kT)
    print("Data set: "+force)
    print("Unweighted mean: {} | std: {}".format(np.mean(data),np.std(data)))
    print("Weighted average: {} | std: {}".format(np.average(data,weights=weights),np.sqrt(np.cov(data,aweights=weights))))   
    print("========================================")