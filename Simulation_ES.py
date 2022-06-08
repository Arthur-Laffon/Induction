# simulation du phénomène d'induction de Faraday

from math import dist, pi, sqrt, atan, acos, asin, sin, cos
import numpy as np

class magnet () : 
    def __init__(points, remanence, velocity) :
        self.points = points
        self.remanence = remanence
        self.velocity = velocity

class coil() : 
    def __init__(position, shape, height, turns) : 
        self.position = position
        self.shape = shape
        self.height = height
        self.N = turns

def faces (points) : 
    #faces = ((bottom/up), (left/right), (front/back))
    faces = [((points[0], points[1], points[2], points[3]), (points[4], points[5], points[6], points[7])), ((points[0], points[1], points[4], points[5]), (points[2], points[3], points[6], points[7])), ((points[0], points[3], points[4], points[7]), (points[1], points[2], points[5], points[6]))]
    return faces

def droites (face) : 
    d1 = (np.subtract(face[0][0], face[1][0]), face[0][0])
    d2 = np.subtract((face[0][1], face[1][1], face[0][1]))
    d3 = np.subtract((face[0][2], face[1][2], face[0][2]))
    d4 = np.subtract((face[0][3], face[1][3], face[0][3]))
    return (d1, d2, d3, d4)

def p_droites (droites, eq) : 
    p = []
    for droite in droites : 
        #define left-hand side of equation
        left_side = np.array(droite, [eq[0], eq[1], eq[2]], np.add(droite, [eq[0], eq[1], eq[2]]))

        #define right-hand side of equation
        right_side = np.array([0, eq[3], eq[3]])

        #solve for x, y, and z
        p.append(np.linalg.inv(left_side).dot(right_side))
    return p

def equation (point, theta, phi) : 
    Cx = point[0]+2*cos(theta)*cos(phi)
    Cy = point[1]+2*sin(theta)*cos(phi)
    Cz = point[2]+2*sin(-phi)
    C = (Cx, Cy, Cz)
    N = np.subtract(C, point)
    d = sum(np.multiply(N, point))
    return (N[0], N[1], N[2], d)

def surface (points) : 
    return dist(points[0],points[1])*dist(points[0],points[3])

def in_plane (point, plane) : 
    x = [plane[i][0] for i in range (4)]
    y = [plane[i][1] for i in range (4)]
    z = [plane[i][2] for i in range (4)]
    if min(x)<point[0]<max(x) and min(y)<point[1]<max(y) and min(z)<point[2]<max(z) : return True
    else : return False

def rect2sphr (vecteur) : 
    """convertit vecteur rectangulaire en coordonnées sphériques"""
    r = sqrt(vecteur[0]**2+vecteur[1]**2+vecteur[2]**2)
    if r==0 : return [0, 0, 0]
    phi = acos(vecteur[2]/r)
    print(sin(phi))
    theta = acos(vecteur[0]/r/sin(phi))
    return [r, theta, phi]

def sphr2rect (vecteur) : 
    """convertit un vecteur sphérique en coordonnées rectangulaires"""
    x = vecteur[0]*sin(vecteur[2])*cos(vecteur[1])
    y = vecteur[0]*sin(vecteur[2])*sin(vecteur[1])
    z = vecteur[0]*cos(vecteur[2])
    return (x, y, z)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def closest (pavé, point) : 
    """calcule le point à la surface du pavé le plus proche de point"""
    faces = faces(pavé)
    spherique = rect2sphr(np.subtract(pavé[3], pavé[0]))
    #sur l'axe x (0 n'a que la coordonnée x)
    theta = spherique[1]
    #sur l'axe z (pi/2 n'a que la coordonnée z)
    phi = spherique[2]
    eq = equation(point, theta, phi)
    for face in faces : 
        p_lignes = p_droites(droites(face), eq)
        d0 = np.subtract(face[0][0], p_lignes[0])
        d1 = np.subtract(face[1][0], p_lignes[0])
        d = min(d0, d1)
        if in_plane(point, p_lignes) : 
            return (np.add(d, point), face, d, theta, phi)

def B (point, magnet) : 
    C = closest(magnet.points, point)
    A = dist(C[1][0][0], C[1][0][1])
    B = dist(C[1][0][0], C[1][0][3])
    L = dist(C[1][0][0], C[1][1][0])
    X = C[1]
    aire = A*B
    magnitude = magnet.remanence/pi*(atan(aire/(2*X*sqrt(4*X**2+A**2+B**2)))-atan(aire/(2*(X+L)*sqrt(4*(X+L)**2+A**2+B**2))))
    v = sphr2rect(magnitude, C[4], C[3])
    return v

def magnetic_flux (magnet, area, resolution) : 
    """evaluate the magnetic flux at a given area"""
    longueur = dist(area[0], area[1])
    largeur = dist(area[0], area[3])
    v_longueur = np.subtract(area[1], area[0])/resolution
    v_largeur = np.subtract(area[3], area[0])/resolution
    surface = surface(area)
    spherique = rect2sphr(np.subtract(area[3], area[0]))
    #sur l'axe x (0 n'a que la coordonnée x)
    theta = spherique[1]
    #sur l'axe z (pi/2 n'a que la coordonnée z)
    phi = spherique[2]
    area_vector = sphr2rect((1, theta, phi))
    points = []
    p = [area[0]]
    B_points = []
    for i in np.arange(0, longueur, longueur/resolution):
        for j in np.arange(0, largeur, largeur/resolution) :
            points.append(p)
            p = np.add(p, v_largeur)
        p = np.sub(p, v_largeur*resolution)
        p = np.add(p, v_longueur)
    for point in points : 
        mag = dist(B(point, magnet), (0, 0, 0))
        rad = angle_between(area_vector, B(point, magnet))
        B_points.append(cos(rad)*mag)
    return sum(B_points)/surface**2