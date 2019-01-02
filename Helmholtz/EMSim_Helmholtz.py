
# coding: utf-8

# In[6]:

x = [1,2,3]


# ## Imports

# In[9]:

import math
import numpy as np
from numba import njit
import numpy.linalg as dLA
import scipy.sparse.linalg as sLA
import scipy.sparse as sp
import scipy.sparse.linalg as sLA
from scipy.sparse.linalg import dsolve
from scipy.sparse import csc_matrix
import itertools
import PIL
from colorize import *


# In[10]:

from scipy.sparse.linalg import dsolve


# ## Notes

# Notes on coodinate systems:
#  - `(x, y)`: For display and specification.  Origin is at 0,0 in the lower left hand corner.
#  - `(i, j)`: Location in matrix. i is row, j is column.  (i, j) = (y, x)

# ## Function Library

# In[11]:

class EMSim:

    def __init__(self, xyShape=(10, 10), WL0=20,
                 epsStart=1. + 0.j, epsEnd=2. + 0.01j,
                 varBox=((0, 0), (0, 0))):
        xMax, yMax = xyShape
        self.xyShape = xyShape
        self.epsStart = epsStart
        self.epsEnd = epsEnd
        self.varBox = varBox
        self.WL0 = WL0  # Size of wl0 in 'pixels'
        self.margin = int(math.ceil(self.WL0))  # Size of PML in 'pixels'
        self.xyFullShape = tuple(
            x + 2 * self.margin for x in self.xyShape)  # shape with PML
        self.k0 = 2 * np.pi / WL0  # k0, free-space avenumber.
        self.initGeoFields()

    def initGeoFields(self):
        m = self.margin
        (xMax, yMax) = self.xyFullShape
        ijShape = (yMax, xMax)

        # The three vars below are all views of the same data block
        # (n,m,3) fields into PML
        self.fullFields = np.zeros(ijShape, np.complex)
        self.fields = self.fullFields[m:-m, m:-m]
        self.fullFields1D = self.fullFields.ravel()  # fields with PML as 1D vector

        # (n,m,3) fields into PML
        self.fullGoal = np.zeros(ijShape, np.complex)
        self.goal = self.fullGoal[m:-m, m:-m]
        self.fullGoalMask = np.zeros(ijShape, np.uint8)
        self.goalMask = self.fullGoalMask[m:-m, m:-m]

        # (n,m,3) fields into PML
        self.fullSources = np.zeros(ijShape, np.complex)
        self.sources = self.fullSources[m:-m, m:-m]
        self.fullSources1D = self.fullSources.ravel()  # fields with PML as 1D vector

        # The three vars below are all views of the same data block
        self.fullEps = np.full(ijShape, self.epsStart,
                               np.complex)  # permittivity matrix
        self.eps = self.fullEps[m:-m, m:-m]  # permittivity matrix
        self.fullEps1D = self.fullEps.ravel()  # permittivity as 1D vector

        ((xmin, xmax), (ymin, ymax)) = self.varBox
        self.varEps = self.fullEps[m + ymin: m + ymax, m + xmin: m + xmax]

        # The three vars below are all views of the same data block
        # permeability matrix
        self.fullMu = np.full(ijShape,  (1 - 0.j), np.complex)
        self.mu = self.fullMu[m:-m, m:-m]
        self.fullMu1D = self.fullMu.ravel()  # permeability as 1D vector

        (self.fullSx, self.fullSy) = self.buildSArrays()
        self.fullSx1D = self.fullSx.ravel()  # permeability as 1D vector
        self.fullSy1D = self.fullSy.ravel()  # permeability as 1D vector

    def buildSArrays(self):
        margin = self.margin
        WL0 = self.WL0

        def s(x):
            return(1 - (1.1) * (1 / WL0) * 6.413j * (x / margin)**4)
        (xDim, yDim) = self.xyShape
        sEdge = np.array([s(x) for x in range(1, margin + 1)])
        middleEdgeX = np.ones(xDim, np.complex)
        middleEdgeY = np.ones(yDim, np.complex)
        ySizeFull = yDim + 2 * margin
        xSizeFull = xDim + 2 * margin
        sEdgeX = np.concatenate((sEdge[::-1], middleEdgeX, sEdge))
        sEdgeY = np.concatenate((sEdge[::-1], middleEdgeY, sEdge))
        sX = np.tile(sEdgeX, (ySizeFull, 1))
        sY = np.tile(sEdgeY, (xSizeFull, 1)).transpose()
        return (sX, sY)

    def setPointSource(self, xy, val):
        (x, y) = xy
        self.sources[y, x] = val

    def setGoalPoint(self, xy, val):
        (x, y) = xy
        self.goal[(y, x)] = val
        self.goalMask[(y, x)] = 1

    def getScore(self):
        score = np.sum(np.abs(self.fields - self.goal)**2 * self.goalMask)
        count = np.sum(self.goalMask)
        aveError = score / count
        return aveError

    def getScoreAbs(self):
        score = np.sum((np.abs(self.fields) - np.abs(self.goal))
                       ** 2 * self.goalMask)
        count = np.sum(self.goalMask)
        aveError = score / count
        return aveError

    def setEpsPoint(self, xy, eps):
        (x, y) = xy
        self.eps[y, x] = eps

    def setMaterialPoint(self, xy, mat):
        (x, y) = xy
        self.eps[y, x] = (1 - mat) * self.epsStart + mat * self.epsEnd

    def setVarMat1D(self, matArray):
        newEps = (1 - matArray) * self.epsStart + matArray * self.epsEnd
        newEps2D = newEps.reshape(self.varEps.shape)
        self.varEps[:] = newEps2D

    def getVarMat1D(self):
        arrayReal = np.real(self.varEps)
        minRange = np.real(self.epsStart)
        maxRange = np.real(self.epsEnd)
        rescaled = (arrayReal - minRange) /             (maxRange - minRange)  # range is now [0, 1]
        return rescaled.ravel()

    def setVarMat2D(self, matArray):
        newEps2D = (1 - matArray) * self.epsStart + matArray * self.epsEnd
        self.varEps[:] = newEps2D

    def getVarMat2D(self):
        arrayReal = np.real(self.varEps)
        minRange = np.real(self.epsStart)
        maxRange = np.real(self.epsEnd)
        rescaled = (arrayReal - minRange) /             (maxRange - minRange)  # range is now [0, 1]
        return rescaled

    def visualizeMaterial(self, maxVal=2, zoom=3):
        minVal = np.real(self.epsStart)
        maxVal = np.real(self.epsEnd)
        pixArray = colorizeArray(self.eps[::-1], min_max=(minVal, maxVal),
                                 colors=([210, 210, 180], [255, 199, 0]),
                                 preFunc=lambda x: np.real(x))
        zoom = 3
        image = PIL.Image.fromarray(pixArray)
        (ySize, xSize) = pixArray.shape[0:2]
        imageBig = image.resize(
            (ySize * zoom, xSize * zoom), PIL.Image.NEAREST)
        return imageBig

    def visualizeSources(self, zoom=3):
        maxRange = np.max(np.abs(self.sources))
        pixArray = colorizeComplexArray(
            self.sources[::-1], maxRad=maxRange, centerColor='black')
        image = PIL.Image.fromarray(pixArray)
        (ySize, xSize) = pixArray.shape[0:2]
        imageBig = image.resize(
            (ySize * zoom, xSize * zoom), PIL.Image.NEAREST)
        return imageBig

    def visualizeGoal(self, zoom=3):
        maxRange = np.max(np.abs(self.goal))
        pixArray = colorizeComplexArray(
            self.goal[::-1], maxRad=maxRange, centerColor='black')
        image = PIL.Image.fromarray(pixArray)
        (ySize, xSize) = pixArray.shape[0:2]
        imageBig = image.resize(
            (ySize * zoom, xSize * zoom), PIL.Image.NEAREST)
        return imageBig

    def visualizeFields(self, zoom=3, maxRange=1):
        maxRange = np.max(np.abs(self.fields))
        pixArray = colorizeComplexArray(
            self.fields[::-1], maxRad=maxRange, centerColor='black')
        image = PIL.Image.fromarray(pixArray)
        (ySize, xSize) = pixArray.shape[0:2]
        imageBig = image.resize(
            (ySize * zoom, xSize * zoom), PIL.Image.NEAREST)
        return imageBig

    def visualizeFieldsMag(self, zoom=3, maxRange=1):
        maxRange = np.max(np.abs(self.fields))
        pixArray = colorizeArray(self.fields[::-1], min_max=(0, maxRange),
                                 colors=([0, 0, 0], [255, 0, 0],
                                         [255, 255, 255]),
                                 preFunc=lambda x: np.abs(x))
        image = PIL.Image.fromarray(pixArray)
        (ySize, xSize) = pixArray.shape[0:2]
        imageBig = image.resize(
            (ySize * zoom, xSize * zoom), PIL.Image.NEAREST)
        return imageBig

    def visualizeFieldsMagWithAbsorber(self, zoom=3, maxRange=1):
        maxRange = np.max(np.abs(self.fullFields))
        pixArray = colorizeArray(self.fullFields[::-1], min_max=(0, maxRange),
                                 colors=([0, 0, 0], [255, 0, 0],
                                         [255, 255, 255]),
                                 preFunc=lambda x: np.abs(x))
        image = PIL.Image.fromarray(pixArray)
        (ySize, xSize) = pixArray.shape[0:2]
        imageBig = image.resize(
            (ySize * zoom, xSize * zoom), PIL.Image.NEAREST)
        return imageBig

    def buildEquations(self):
        (xMax, yMax) = self.xyFullShape
        fShape = self.xyFullShape
        k0 = self.k0
        (rList, cList, vList) = allEqGen(fShape, self.fullEps, self.fullMu,
                                         self.fullSx, self.fullSy, k0)
        n = self.fullFields.size
        self.M = sp.coo_matrix((vList, (rList, cList)), (n, n))

    def solve(self):
        B = self.fullSources1D
        X = dsolve.spsolve(csc_matrix(self.M), B)
        self.fullFields1D[:] = X


# In[12]:

@njit
def EzInd(xy, xyFullShape):
    (xMax, yMax) = xyFullShape
    (x, y) = xy
    (xWrap, yWrap) = (x % xMax, y % yMax)
    (i, j) = (yWrap, xWrap)
    return j + (i * xMax)


# In[13]:

@njit
def MakeEzEq(xy, shape, eps, mu, sx, sy, k0):
    """
    Generates an equation for the discretized form of 
        Lap(E) - k0^2 mu eps E = 0
    Note, when eps = PEC, signified by any im(eps) <= -1000, the equation becomes
        E = 0
    """
    eqNum = EzInd(xy, shape)
    (x, y) = xy
    delta = (1. + 0.j)
    # Calculate Coefficients
    cEz = -1. + 0j  # when brought to the RHS of the EQ
    k = (eps**0.5 * mu**0.5) * k0

    cX = 1 / ((k0 * sx * delta)**2 * eps * mu)
    cY = 1 / ((k0 * sy * delta)**2 * eps * mu)

    (cEzp0, cEzn0, cEz0p, cEz0n, cEz00) = (
        cX, cX, cY, cY, -2 * cX - 2 * cY + 1)
    # Calculate Indices
    (iEzp0, iEzn0, iEz0p, iEz0n, iEz00) = (EzInd((x + 1, y), shape), EzInd((x - 1, y), shape),
                                           EzInd(
                                               (x, y + 1), shape), EzInd((x, y - 1), shape),
                                           EzInd((x, y), shape))
    # Build and return sparse equation additions
    rowArray = (eqNum, eqNum, eqNum, eqNum, eqNum)
    colArray = (iEzp0, iEzn0, iEz0p, iEz0n, iEz00)
    coeffArray = (cEzp0 / sy, cEzn0 / sy, cEz0p / sx, cEz0n / sx, cEz00)
    return (rowArray, colArray, coeffArray)


# In[14]:

@njit
def allEqGen(fShape, epsArr, muArr, sxArr, syArr, k0):
    rowCoeffs = []
    colCoeffs = []
    valCoeffs = []
    (xMax, yMax) = fShape
    for y in range(yMax):
        for x in range(xMax):
            xy = (x, y)
            eps = epsArr[y, x]
            mu = muArr[y, x]
            sx = sxArr[y, x]
            sy = syArr[y, x]
            (newRow0, newCol0, newVal0) = MakeEzEq(
                xy, fShape, eps, mu, sx, sy, k0)
            rowCoeffs.extend(newRow0)
            colCoeffs.extend(newCol0)
            valCoeffs.extend(newVal0)
    return (rowCoeffs, colCoeffs, valCoeffs)


# 

# In[ ]:



