from scipy.constants import mu_0, hbar
import spindata
import numpy as np

# hyperfine Eq. 3.19 in Wenckebach, truncated
def A(index, nuc, r, theta, phi=0):
    gS = spindata.gamma("E")
    gI = spindata.gamma(nuc)
    
    A0 = mu_0 / (4*np.pi)*(hbar*gS*gI)/(r**3)
    
    myDict = {"0" : A0, 
              "z+" : -3*A0*np.cos(theta)*np.sin(theta)*np.exp(1j*phi),
              "++" : -3*A0*(np.sin(theta)**2*np.exp(1j*2*phi)),
              "zz" : A0*(1 - 3*np.cos(theta)**2)}
    
    myDict["z-"] = np.conj(myDict["z+"])
    
    return myDict[index]
