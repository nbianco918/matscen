import numpy as np

h = 6.626075e-34 #J-s
e = 1.602177e-19 #C
m_e = 9.109389e-31 #kg
m_n = 1.674929e-27 #kg
c = 2.99792458e8 #m/s
kb = 1.38e-23 #J/K

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

def metric_matrix(a, b, c, alpha, beta, gamma):
    """
    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.

    Returns
    -------
    g : TYPE
        DESCRIPTION.

    """
    g = np.array([[a**2, a*b*np.cos(np.deg2rad(gamma)), a*c*np.cos(np.deg2rad(beta))],
                  [b*a*np.cos(np.deg2rad(gamma)), b**2, b*c*np.cos(np.deg2rad(alpha))],
                  [c*a*np.cos(np.deg2rad(beta)),c*b*np.cos(np.deg2rad(alpha)) ,c**2]])
    return g

def crystal_dot_product(v1, v2, a, b, c, alpha, beta, gamma):
    """
    Parameters
    ----------
    v1 : TYPE
        DESCRIPTION.
    v2 : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    g = metric_matrix(a, b, c, alpha, beta, gamma)
    return np.dot(v1, np.dot(g, v2.T))

def crystal_length(v1, v2, a, b, c, alpha, beta, gamma):
    """
    Parameters
    ----------
    v1 : TYPE
        DESCRIPTION.
    v2 : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    return np.sqrt(crystal_dot_product(v1, v2, a, b, c, alpha, beta, gamma))

def crystal_angle(v1, v2, a, b, c, alpha, beta, gamma):
    """
    Parameters
    ----------
    v1 : TYPE
        DESCRIPTION.
    v2 : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    
    cos_theta = crystal_dot_product(v1, v2, a, b, c, alpha, beta, gamma)/(crystal_length(v1, v1, a, b, c, alpha, beta, gamma)*crystal_length(v2, v2, a, b, c, alpha, beta, gamma))
    return np.rad2deg(np.arccos(cos_theta))

def volume_of_cell(a, b, c):
    """
    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    Returns
    -------
    V : TYPE
        DESCRIPTION.
    """
    V = np.dot(a, np.cross(b, c))
    return V

def F(alpha, beta, gamma):
    F = np.cos(np.deg2rad(alpha))*np.cos(np.deg2rad(beta))-np.cos(np.deg2rad(gamma))
    return F
def reciprocal_metric_matrix(a, b, c, alpha, beta, gamma, cs):
    if cs == 'cubic':
        g_star = np.array([[1/a**2, 0, 0],
                           [0, 1/a**2, 0],
                           [0, 0, 1/a**2]])
    elif cs == "tetragonal":
        g_star = np.array([[1/a**2, 0, 0],
                           [0, 1/a**2, 0],
                           [0, 0, 1/c**2]])
    elif cs == 'orthorhombic':
        g_star = np.array([[1/a**2, 0, 0],
                           [0, 1/b**2, 0],
                           [0, 0, 1/b**2]])
    elif cs == 'hexagonal':
        g_star = np.array([[4/(3*a**2), 2/(3*a**2), 0],
                           [2/(3*a**2), 4/(3*a**2), 0],
                           [0, 0, 1/c**2]])
    elif cs == 'rhombohedral':
        W_square = a**2*(1+np.cos(np.deg2rad(alpha))**2-2*np.cos(np.deg2rad(alpha))**2)
        g_star = (1/W_square)*np.array([[1+np.cos(np.deg2rad(alpha)), -1*np.cos(np.deg2rad(alpha)), -1*np.cos(np.deg2rad(alpha))],
                           [-1*np.cos(np.deg2rad(alpha)), 1+np.cos(np.deg2rad(alpha)), -1*np.cos(np.deg2rad(alpha))],
                           [-1*np.cos(np.deg2rad(alpha)), -1*np.cos(np.deg2rad(alpha)), 1+np.cos(np.deg2rad(alpha))]])
    elif cs == 'monoclinic':
        g_star = np.array([[1/(a**2*np.sin(np.deg2rad(beta))**2), 0, -1*np.cos(np.deg2rad(beta))/(a*c*np.sin(np.deg2rad(beta))**2)],
                           [0, 1/b**2, 0],
                           [-1*np.cos(np.deg2rad(beta))/(a*c*np.sin(np.deg2rad(beta))**2), 0, 1/(c**2*np.sin(np.deg2rad(beta))**2)]])
    elif cs == 'triclinic':
        V_square = a**2*b**2*c**2*(1-np.cos(np.deg2rad(alpha))**2-np.cos(np.deg2rad(beta))**2-np.cos(np.deg2rad(gamma))**2+2*np.cos(np.deg2rad(alpha))*np.cos(np.deg2rad(beta))*np.cos(np.deg2rad(gamma)))
        
        g_star = (1/V_square)*np.array([[b**2*c**2*np.sin(np.deg2rad(alpha))**2, a*b*c**2*F(alpha, beta, gamma), a*b**2*c*F(gamma, alpha, beta)],
                           [a*b*c**2*F(alpha, beta, gamma), a**2*c**2*np.sin(np.deg2rad(beta))**2, a**2*b*c*F(beta, gamma, alpha)],
                           [a*b**2*c*F(gamma, alpha, beta), a**2*b*c*F(beta, gamma, alpha), a**2*b**2*np.sin(np.deg2rad(gamma))**2]])
    return g_star

def reciprocal_crystal_dot_product(v1, v2, a, b, c, alpha, beta, gamma, cs):
    g_star = reciprocal_metric_matrix(a, b, c, alpha, beta, gamma, cs)
    return np.dot(v1, np.dot(g_star, v2.T))

def reciprocal_crystal_length(v1, v2, a, b, c, alpha, beta, gamma, cs):
    
    return np.sqrt(reciprocal_crystal_dot_product(v1, v2, a, b, c, alpha, beta, gamma, cs))

def interplanar_spacing(v1, v2, a, b, c, alpha, beta, gamma, cs):
    g_mag = reciprocal_crystal_length(v1, v2, a, b, c, alpha, beta, gamma, cs)

    return 1/g_mag

def reciprocal_crystal_angle(v1, v2, a, b, c, alpha, beta, gamma, cs): 
    
    cos_theta = reciprocal_crystal_dot_product(v1, v2, a, b, c, alpha, beta, gamma, cs)/(reciprocal_crystal_length(v1, v1, a, b, c, alpha, beta, gamma, cs)*reciprocal_crystal_length(v2, v2, a, b, c, alpha, beta, gamma, cs))
    return np.rad2deg(np.arccos(cos_theta))

d = interplanar_spacing(np.array([1,1,0]), np.array([1,1,0]), 2, 2, 2, 90, 90, 90, 'cubic')
print(d)

d = interplanar_spacing(np.array([1,1,1]), np.array([1,1,1]), 3, 4, 6, 90, 90, 120, 'triclinic')
print(d) 

print(reciprocal_crystal_angle(np.array([1,0,1]), np.array([-2,0,1]), 4, 6, 5, 90, 120, 90, 'monoclinic'))

