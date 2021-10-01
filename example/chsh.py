from julia import SPDC
import numpy as np
from scipy.optimize import minimize

def chsh(p,η=1.0):
    """ Compute the CHSH value for parameter p

    Param
    -----
        p : np.array 
            Circuit parameters, 1-d float array.
        η : float
            Photondetector efficiency

    Return
    ------
       chsh_score : float
            The computed CHSH score
    """
    X = 2 # Number of measurement of Alice
    Y = 2 # Number of measurement of Bob
    # Translate the parameter to a NamedTuple
    p = SPDC.param(p,X,Y,ηA=η,ηB=η)
    # Compute all the correlators <A_xB_y>
    AB = SPDC.spdc_correlators(p,X,Y)
    # Compute the CHSH score
    chsh_score = 0.0
    for x in range(2):
        for y in range(2):
            chsh_score += (-1)**(x*y)*AB[x,y] 
    return chsh_score

def optimize_chsh(η=1.0,x0=np.random.random(10)):
    """ Optimize the CHSH score. """
    res = minimize(lambda x: -chsh(x,η=η),
            x0,
            method="BFGS")
    return res

def chsh_with_loss(step=1e-2):
    """ Optimized CHSH score for decreasing efficiency value """
    scores = []
    η = 1.0
    r = optimize_chsh()
    while -r.fun >= 2.0:
        scores.append([η,r.fun])
        η -= step
        r = optimize_chsh(η=η,x0=r.x)
        print(scores[-1])
    return np.array(scores)
