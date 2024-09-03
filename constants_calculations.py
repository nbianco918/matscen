import numpy as np

# Define constants
h = 6.626075e-34 #J-s - Planck's constant
e = 1.602177e-19 #C - charge of an electron
m_e = 9.109389e-31 #kg - mass of an electron
m_n = 1.674929e-27 #kg - mass of a neutron
c = 2.99792458e8 #m/s - speed of light
kb = 1.38e-23 #J/K - Boltzmann's constat

def energy_photon(nu = "", lambda_ = ""):
    """
    Parameters
    ----------
    nu : TYPE, optional
        DESCRIPTION. The default is "".
    lambda_ : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    
    if nu == "":
        return ((h*c)/lambda_)/e
    else:
        return (h*nu)/e
    
def calculate_frequency(lambda_):
    """
    Parameters
    ----------
    lambda_ : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    """
    c = lambda * nu 
    """
    return c/lambda_

def wavelength_accelerated_electron(V):
    """
    Parameters
    ----------
    V : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    return (h/np.sqrt(2*m_e*e*V))

def wavelength_accelerated_electron_rc(V):
    """
    Parameters
    ----------
    V : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    return (h/np.sqrt(2*m_e*e*V*(1+(e*V)/(2*m_e*c**2))))
    
def deBroglie_wavelength(m, v):
    """
    Parameters
    ----------
    m : TYPE
        DESCRIPTION.
    v : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return h/(m*v)

def deBroglie_wavelength_thermal_neutrons(m, T):
    """
    Parameters
    ----------
    m : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    return h/(np.sqrt(3*m*kb*T))

def temp_of_n_given_lambda(lambda_):
    return (h**2)/(3*m_n*kb*lambda_**2)

def hydrogenic_atom(Z, n):
    """
    Parameters
    ----------
    Z : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    # returns eV
    return 13.606*(Z**2/n**2)

def atomic_radii_bcc(a):
    """
    Parameters
    ----------
    a : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    return (np.sqrt(3)*a)/4

def atomic_radii_equiatomic_bcc(a, r1):
    """
    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    r1 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    return (np.sqrt(3)*a - 2*r1)/2

def fcc_atom_density(a, V_mat):
    """
    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    V_mat : TYPE
        DESCRIPTION.

    Returns
    -------
    N_total : TYPE
        DESCRIPTION.
    """
    # 2 atoms per cell
    cell_density = 2/(a**3) # atoms/cell volume (m3)
    N_total = V_mat*cell_density
    return N_total