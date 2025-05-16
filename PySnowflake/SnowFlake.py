######################################################################
#
#   Snowflake2 for python3
#           A.Vladimirov (16.05.2025)
#
######################################################################

import pysnowflake
import numpy

def initialize(fileName):
    """
    Initialization of snowflake

    Parameters
    ----------
    fileName : string
        The path to the constants-file

    Returns
    -------
    None.

    """
    
    import os.path
    if os.path.exists(fileName):
    
        pysnowflake.pysnowflake.initialize(fileName)
    else:
        raise FileNotFoundError('INI-file '+fileName+' NOT FOUND')

def UpdateTables(mu0,mu1):
    """
    

    Parameters
    ----------
    mu0 : float
        Initial scale of the evolution tables. At this scale the boundary is defined.
    mu1 : float
        Max value of scale

    Returns
    -------
    None.

    """
    
    if isinstance(mu0,float) and isinstance(mu1,float):
        if(mu0<mu1):
            pysnowflake.pysnowflake.updateevolutiontable(mu0,mu1)
        else:
            raise ValueError("mu0 should be smaller than mu1")
    else:
        raise TypeError("mu0 and mu1 Must be float")
            
def setNPparameters(l):
    """
    Set NP parameters 

    Parameters
    ----------
    l : float array
        list of NP parameters

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        pysnowflake.pysnowflake.updatenpparameters(numpy.asfortranarray(l))
    else:
        raise TypeError()
        
def G2List(x,Q,f):
    return pysnowflake.pysnowflake.snowflake_g2_list(numpy.asfortranarray(x),\
                                       numpy.asfortranarray(Q),\
                                       numpy.asfortranarray(f),
                                       len(x))
        
def D2List(Q,f):
    return pysnowflake.pysnowflake.snowflake_d2_list(numpy.asfortranarray(Q),\
                                       numpy.asfortranarray(f),
                                       len(Q))