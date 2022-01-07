from math import *
import scipy
from sympy import *
import scipy.integrate

a = 6378249.145
b = 6356515
f = 1 - (b / a)


def exc():
    return 1 - (a ** 2 / b ** 2)


def exc1(a, b):
    return 1 - (a ** 2 / b ** 2)


# Fonction 1 Parametres d'un ellipse

def param_1(a, b):
    d = dict()
    d['f'] = round(float((1 - b / a)), 8)
    d['e1'] = round(float(1 - (pow(b, 2) / pow(a, 2))), 8)
    d['e2'] = round(float((pow(a, 2) / pow(b, 2)) - 1), 8)
    d['alpha'] = round(float(degrees(acos(b / a))), 8)
    d['c'] = round(float(pow(a, 2) / b), 8)
    return d


# Fonction 2

def param_phi(X, Y, Z, a, b):
    phi = dict()
    R = sqrt((X ** 2) + (Y ** 2) + (Z ** 2))
    r = sqrt((X ** 2) + (Y ** 2))
    u = atan((Z / r) * ((1 - f) + ((exc1(a, b) * a) / R)))
    phi_rad = atan((Z * (1 - f) + exc1(a, b) * a * (sin(u)) ** 3) / ((1 - f) * (r - exc1(a, b) * a * (cos(u) ** 3))))
    phi_deg = degrees(phi_rad)
    phi1 = floor(phi_deg)
    phi2 = floor((phi_deg - phi1) * 60)
    phi3 = floor((((phi_deg - phi1) * 60) - phi2) * 60)
    phi1 = round(phi1, 2)
    phi2 = round(phi2, 2)
    phi3 = round(phi3, 2)
    phi['p1'] = str(phi1)
    phi['p2'] = str(phi2)
    phi['p3'] = str(phi3)

    return phi


def param_lamda(X, Y):
    d = dict()

    lamda = degrees(atan(Y / X))
    lamda1 = floor(lamda)
    lamda2 = floor((lamda - lamda1) * 60)
    lamda3 = floor((((lamda - lamda1) * 60) - lamda2) * 60)
    lamda1 = round(lamda1, 2)
    lamda2 = round(lamda2, 2)
    lamda3 = round(lamda3, 2)
    d['l1'] = str(lamda1)
    d['l2'] = str(lamda2)
    d['l3'] = str(lamda3)

    return d


def param_h(X, Y, Z, a, b):
    R = sqrt((X ** 2) + (Y ** 2) + (Z ** 2))
    r = sqrt((X ** 2) + (Y ** 2))
    u = atan((Z / r) * ((1 - f) + ((exc1(a, b) * a) / R)))
    phi_rad = atan((Z * (1 - f) + exc1(a, b) * a * (sin(u)) ** 3) / ((1 - f) * (r - exc1(a, b) * a * (cos(u) ** 3))))
    h = (r * cos(phi_rad)) + (Z * sin(phi_rad)) - (a * sqrt(1 - (exc1(a, b)) * (sin(phi_rad) ** 2)))
    h = round(Abs(h), 2)
    return h


# Calcul de rayon de courbure du premier vertical
def Rayon_vertical(phi, a, b):
    def get_w(phi_w):
        w = sqrt(1 - exc1(a, b) * (sin(phi_w) ** 2))
        return w

    phi_rad = radians(phi)
    W = get_w(phi_rad)
    N = a / W
    N = round(N, 2)
    return N


# coor géographiques to coor cartèsiennes
def calcul_X(phi1, phi2, phi3, lamda1, lamda2, lamda3, h, a, b):
    phi = phi1 + phi2 / 60 + phi3 / 3600
    lamda = lamda1 + lamda2 / 60 + lamda3 / 3600
    lamda_rad = radians(lamda)
    phi_rad = radians(phi)
    x = ((Rayon_vertical(phi, a, b) + h) * cos(phi_rad) * cos(lamda_rad))
    x = round(x, 4)
    return x


def calcul_Z(phi1, phi2, phi3, h, a, b):
    phi = phi1 + phi2 / 60 + phi3 / 3600
    phi_rad = radians(phi)
    Z = ((Rayon_vertical(phi, a, b) * (1 - exc1(a, b))) + h) * sin(phi_rad)
    Z = round(Z, 4)
    return Z


def calcul_Y(phi1, phi2, phi3, lamda1, lamda2, lamda3, h, a, b):
    phi = phi1 + phi2 / 60 + phi3 / 3600
    lamda = lamda1 + lamda2 / 60 + lamda3 / 3600
    lamda_rad = radians(lamda)
    phi_rad = radians(phi)
    Y = (Rayon_vertical(phi, a, b) + h) * cos(phi_rad) * sin(lamda_rad)
    Y = round(Y, 4)
    return Y


# Fonction 3

# Calcul de rayon de courbure du meridien
def Rayon_meridien(phi):
    def get_w(phi_w):
        w = sqrt(1 - exc() * sin(phi_w))
        return w

    phi_rad = radians(phi)
    W = get_w(phi_rad)
    M = a * (1 - exc()) / W ** 3
    M = round(float(M), 2)
    return M


# Calcul de rayon de courbure du premier vertical
def Rayon_vertical2(phi):
    def get_w(phi_w):
        w = sqrt(1 - exc() * sin(phi_w))
        return w

    phi_rad = radians(phi)
    W = get_w(phi_rad)
    N = a / W
    N = round(float(N), 2)
    return N


# Calcul de rayon de courbure de section normal d'azimut
def Rayon_azimut(alpha, phi):
    phi_rad = radians(phi)
    R = (Rayon_meridien(phi_rad) * sin(alpha) ** 2) + (Rayon_vertical2(phi_rad) * cos(alpha) ** 2) / Rayon_vertical2(
        phi_rad) * Rayon_meridien(phi_rad)
    R = round(float(R), 2)
    return R


# Fonction 4
# phi
def calcul_beta_phi(phi):
    phi_rad = radians(phi)
    beta = float(atan(b / a * tan(phi_rad)))
    beta_deg = degrees(beta)
    beta_deg = round(float(beta_deg), 4)
    return beta_deg


def calcul_psy_phi(phi):
    phi_rad = radians(phi)
    psy = float(atan(b ** 2 / a ** 2 * tan(phi_rad)))
    psy_deg = degrees(psy)
    psy_deg = round(float(psy_deg), 4)
    return psy_deg


# psy
def calcul_phi_psy(psy):
    psy_rad = radians(psy)
    phi = float(atan(a / b * tan(psy_rad)))
    phi_deg = degrees(phi)
    phi_deg = round(float(phi_deg), 4)
    return phi_deg


def calcul_beta_psy(psy):
    psy_rad = radians(psy)
    beta = float(atan(a ** 2 / b ** 2 * tan(psy_rad)))
    beta_deg = degrees(beta)
    beta_deg = round(float(beta_deg), 4)
    return beta_deg


# beta

def calcul_psy_beta(beta):
    beta_rad = radians(beta)
    psy = float(atan(b / a * tan(beta_rad)))
    psy_deg = degrees(psy)
    psy_deg = round(float(psy_deg), 4)
    return psy_deg


def calcul_phi_beta(beta):
    beta_rad = radians(beta)
    phi = float(atan(a / b * tan(beta_rad)))
    phi_deg = degrees(phi)
    phi_deg = round(float(phi_deg), 4)
    return phi_deg


# Fonction 5
# CALCUL DE LONGUEUR D'UN ARC DE MERIDIEN
def calcul_Arc_meridien(phi1, phi2):
    # CALCUL DE A,B,C,D,E,F
    A = 1 + (3 / 4) * exc() ** 2 + (45 / 64) * exc() ** 4 + (175 / 256) * exc() ** 6 + (11025 / 16348) * exc() ** 8 \
        + (43659 / 65536) * exc() ** 10
    B = (3 / 4) * exc() ** 2 + (15 / 16) * exc() ** 4 + (525 / 512) * exc() ** 6 + (2205 / 2048) * exc() ** 8 + (
            72765 / 65536) * exc() ** 10
    C = (15 / 16) * exc() ** 4 + (105 / 256) * exc() ** 6 + (2205 / 4096) * exc() ** 8 + (10395 / 16348) * exc() ** 10
    D = (35 / 512) * exc() ** 6 + (315 / 2048) * exc() ** 8 + (10185 / 131072) * exc() ** 10
    E = (315 / 16384) * exc() ** 8 + (3465 / 85336) * exc() ** 10
    F = (639 / 131072) * exc() ** 10

    # calcul des parametres
    Z1 = A * a * (1 - exc())
    Z2 = (B / 2) * a * (1 - exc())
    Z3 = (C / 4) * a * (1 - exc())
    Z4 = (D / 6) * a * (1 - exc())
    Z5 = (E / 8) * a * (1 - exc())
    Z6 = (F / 10) * a * (1 - exc())

    phi1C = radians(phi1)
    phi2C = radians(phi2)

    M1 = Z1 * phi1C - Z2 * sin(2 * phi1C) + Z3 * sin(4 * phi1C) - Z4 * sin(6 * phi1C) + Z5 * sin(
        8 * phi1C) - Z6 * sin(10 * phi1C)
    M2 = Z1 * phi2C - Z2 * sin(2 * phi2C) + Z3 * sin(4 * phi2C) - Z4 * sin(6 * phi2C) + Z5 * sin(
        8 * phi2C) - Z6 * sin(10 * phi2C)

    M = abs(M1 - M2)
    M = round(float(M), 2)
    return M


def calcul_Arc_parallele(lamda1, lamda2, phi):
    phiC = radians(phi)
    w = sqrt(1 - exc() * (pow(sin(phiC), 2)))
    l = a / w
    R = l * cos(phiC)
    N = abs(R * (lamda1 - lamda2))
    N = round(float(N), 2)
    return N


# Fonction 6 Surface partie terreste


def calculeSurface(phi1, phi2, lamda1, lamda2):
    phi1C = radians(phi1)
    phi2C = radians(phi2)
    lamda1C = radians(lamda1)
    lamda2C = radians(lamda2)

    fun = lambda x: cos(x) * (1 - exc() * sin(x) ** 2) ** -2
    integrale = scipy.integrate.quad(fun, phi1C, phi2C)
    Surface = float(b ** 2 * (lamda2C - lamda1C) * integrale[0])
    Surface = round(float(Surface), 2)
    return Surface


# Fonction 7 : PB Direct
# Calcul de sigma12
def calcul_sigma(d):
    R = d / (((2 * a) + b) / 3)
    # sigma12 = d / R
    return R


def calcul_phi2(phi1, d, A12):
    phi1C = radians(phi1)
    A12 = radians(A12)
    phi2 = asin(sin(phi1C) * cos(calcul_sigma(d)) + cos(phi1C) * sin(calcul_sigma(d)) * cos(A12))
    phi2C = degrees(phi2)
    phi2C = round(phi2C, 4)
    return phi2C


def calcul_lamda2(lamda1, phi1, A12, d):
    phi1C = radians(phi1)
    A12 = radians(A12)
    lamda1C = radians(lamda1)
    lamda2 = lamda1C + acot((cot(calcul_sigma(d)) * cos(phi1C) - sin(phi1C) * cos(A12)) / sin(A12))
    lamda2C = degrees(lamda2)
    lamda2C = round(lamda2C, 4)
    return lamda2C


def calcul_A21_direct(phi1, A12, d):
    phi1C = radians(phi1)
    A12C = radians(A12)
    A21 = acot((cos(calcul_sigma(d)) * cos(A12C) - tan(phi1C) * sin(calcul_sigma(d))) / sin(A12C))
    A21C = degrees(A21)
    A21C = round(A21C, 4)
    return A21C


# Fonction 8 : PB INVERSE

def calcul_sigma12(phi1, phi2, lamda1, lamda2):
    phi1C = radians(phi1)
    phi2C = radians(phi2)
    lamda1C = radians(lamda1)
    lamda2C = radians(lamda2)
    sigma12 = acos((sin(phi1C) * sin(phi2C)) + (cos(phi1C) * cos(phi2C) * cos(lamda2C - lamda1C)))
    sigma12 = round(float(sigma12), 6)
    return sigma12


def calcul_A12_inverse(phi1, phi2, lamda1, lamda2):
    phi1C = radians(phi1)
    phi2C = radians(phi2)
    lamda1C = radians(lamda1)
    lamda2C = radians(lamda2)
    A12 = float(acot(((tan(phi2C) * cos(phi1C)) - (sin(phi1C) * cos(lamda2C - lamda1C)) / sin(lamda2C - lamda1C))))
    A12C = degrees(A12)
    A12C = round(float(A12C), 6)
    return A12C


def calcul_A21_inverse(phi1, phi2, lamda1, lamda2):
    phi1C = radians(phi1)
    phi2C = radians(phi2)
    lamda1C = radians(lamda1)
    lamda2C = radians(lamda2)
    A21 = -acot((tan(phi1C) * cos(phi2C)) - (sin(phi2C) * cos(lamda2C - lamda1C)) / sin(lamda2C - lamda1C))
    A21C = degrees(A21)
    A21C = round(float(A21C), 6)
    return A21C
