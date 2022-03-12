import numpy as np
from math import log
import matplotlib.pyplot as plt

class Lyapunov():
  def __init__(self,IOTA = 10e-15, MAX_LN_R = 12, MIN_LN_R = -12, N_LN_R = 600):
    self.IOTA = IOTA
    self.MAX_LN_R = MAX_LN_R
    self.MIN_LN_R = MIN_LN_R
    self.N_LN_R = N_LN_R

  def AllocateDMatrix(self,nRows, nCols):
    return np.zeros((nRows,nCols))

  def ComputeSlopes(self,gMaxDivergeT,gDivergence,gNTests,gCSum):
    IOTA = self.IOTA
    MAX_LN_R = self.MAX_LN_R
    MIN_LN_R = self.MIN_LN_R
    N_LN_R = self.N_LN_R
    m = 0
    b = 0
    rr = 0
    if N_LN_R > gMaxDivergeT:
      siz = N_LN_R
    else:
      siz = gMaxDivergeT
    data = np.zeros((siz))
    for i in range(gNTests):
      i2 = i +gNTests
      k = N_LN_R/(MAX_LN_R-MIN_LN_R)
      for j in range(N_LN_R):
        data[j] = gCSum[j][i]
      for j in range(3,N_LN_R-3):
        m,b,rr = self.LineFit(data+j-3, 7,m,b,rr)
        gDivergence[j][i2] = m
      m,b,rr = self.LineFit(data, 5, m, b, rr)
      gDivergence[2][i2] = m
      m,b,rr = self.LineFit(data+gMaxDivergeT-5, 5, m, b, rr)
      gDivergence[gMaxDivergeT-3][i2] = m
      m,b,rr = self.LineFit(data, 3, m, b, rr)
      gDivergence[1][i2] = m
      m,b,rr = self.LineFit(data+gMaxDivergeT-3, 3, m, b, rr)
      gDivergence[gMaxDivergeT-2][i2] = m
      gDivergence[0][i2] = data[1]-data[0]
      gDivergence[gMaxDivergeT-1][i2] = data[gMaxDivergeT-1]-data[gMaxDivergeT-2]

    return gDivergence

  def FreeDMatrix(self,mat, nRows): 
    del mat

  def GenerateTemplateFile(self):
    # FILE   *outFile;
    # outFile = fopen("sample.l1d2", "w")
    # fclose(outFile);
    pass

  def GetData(self,gData):
    return gData, len(gData)

  def LineFit(self,data, n, m, b, rr):
    IOTA = self.IOTA
    MAX_LN_R = self.MAX_LN_R
    MIN_LN_R = self.MIN_LN_R
    N_LN_R = self.N_LN_R
    sx = sy = sxy = sx2 = sy2 = 0
    for i in range(n):
        x = i
        y = data[i]
        sx += x
        sy += y
        sx2 += x*x
        sy2 += y*y
        sxy += x*y
    k = n*sx2-sx*sx
    mTemp = (n*sxy-sx*sy)/k
    bTemp = (sx2*sy-sx*sxy)/k
    k = sy*sy/n
    if k==sy2:
        rrTemp = 1.0
    else:
        rrTemp = (bTemp*sy+mTemp*sxy-k)/(sy2-k)
        rrTemp = 1.0 - (1.0-rrTemp)*(n-1.0)/(n-2.0)
    m = mTemp
    b = bTemp
    rr = rrTemp

    return m, b, rr

  def PercentDone(self,percentDone):
    last=100
    if percentDone<last:
        last = 0
    elif percentDone>last and percentDone%2==0:
        last = percentDone
        if percentDone%10==0:
          print(percentDone/10)
        else:
          print(".")
        
  def ProcessTest(self,testN,gTest,gMaxDivergeT,gNDivergence,gDivergence,gCSum):
    IOTA = self.IOTA
    MAX_LN_R = self.MAX_LN_R
    MIN_LN_R = self.MIN_LN_R
    N_LN_R = self.N_LN_R
    nPts = 0
    nCompletedPairs=0
    m = gTest[testN]['m']
    J = gTest[testN]['J']
    W = gTest[testN]['W']
    divergeT = gTest[testN]['divergeT']
    gData, nPts = self.GetData(gTest[testN]['seriesN'])
    
    k1 = N_LN_R/(MAX_LN_R-MIN_LN_R)
    k1 *= 0.5
    k2 = N_LN_R/2

    nVectors = nPts-J*(m-1)
    
    isNeighbor = [0]*nVectors
    
    for i in range(gMaxDivergeT):
      gNDivergence[i] = gDivergence[i][testN] = 0

    i = 0
    while(i<nVectors):
        percentDone = 100.0*nCompletedPairs/nVectors*2+0.5
        percentDone = 100.0*i/nVectors+0.5
        self.PercentDone(percentDone)
        
        if isNeighbor[i] == 0:
          distance = 10e10
          if i>W:
            for j in range(i-W):
              d=0
              for k in range(m):
                  T = k*J
                  temp = gData[i+T]-gData[j+T]
                  temp *= temp
                  d += temp
              d += IOTA
              CSumIndex = k1*log(d)+k2
              if CSumIndex<0:
                CSumIndex = 0
              if CSumIndex>=N_LN_R:
                CSumIndex = N_LN_R-1

              gCSum[int(CSumIndex)][testN] += 1

              if d<distance:
                  distance=d
                  neighborIndex=j
          if i<nVectors-W:
            for j in range(i+W,nVectors):
              d=0
              for k in range(m):
                  T = k*J
                  temp = gData[i+T]-gData[j+T]
                  temp *= temp
                  d += temp
              d += IOTA

              CSumIndex = k1*log(d)+k2
              if CSumIndex<0:
                CSumIndex = 0
              if CSumIndex>=N_LN_R:
                CSumIndex = N_LN_R-1

              gCSum[int(CSumIndex)][testN] += 1

              if d<distance:
                  distance=d;
                  neighborIndex=j
          isNeighbor[neighborIndex] = 1
          
          for j in range(divergeT+1):
            maxIndex = nPts-m*J-j-1
            if i<maxIndex and neighborIndex<maxIndex:
              d = 0
              for k in range(m):
                  T = k*J+j
                  temp = gData[i+T]-gData[neighborIndex+T]
                  temp *= temp
                  d += temp
              d += IOTA
              gNDivergence[j] += 1
              temp = 0.5*log(d)
              gDivergence[j][testN] += temp
          nCompletedPairs += 1
        i = i + 1

    for i in range(1,N_LN_R):
      gCSum[i][testN] += gCSum[i-1][testN]

    kNorm = 1.0/gCSum[N_LN_R-1][testN]
    for i in range(N_LN_R):
      gCSum[i][testN] *= kNorm

    for i in range(N_LN_R):
      temp = gCSum[i][testN]
      if (temp<0.000045) or (temp>0.990050):
        gCSum[i][testN] = 0
      else:
        gCSum[i][testN] = log(temp)

    for i in range(divergeT+1):
      if gNDivergence[i]>0:
        gDivergence[i][testN] /= gNDivergence[i]

    del isNeighbor
    del gData

    return gDivergence, gCSum

  def lyap_r(self,gTest,gNTests,gMaxDivergeT,gDivergence,debug_plot=False):
    IOTA = self.IOTA
    MAX_LN_R = self.MAX_LN_R
    MIN_LN_R = self.MIN_LN_R
    N_LN_R = self.N_LN_R
    aa = []
    bb = []
    cc = []
    for i in range(gMaxDivergeT):
        for testN in range(gNTests):
          if i<=gTest[testN]['divergeT']:
            aa.append(gDivergence[i][testN])
        rs = testN
        for testN in range(rs,2*gNTests-1):
          if i<=gTest[testN-gNTests]['divergeT']:
            bb.append(gDivergence[i][testN])
        if i<=gTest[testN-gNTests]['divergeT']:
          cc.append(gDivergence[i][testN])
    
    alpha = np.polyfit(np.arange(gTest[0]['divergeT']+1),aa[0:],1)

    if debug_plot == True:
        plt.figure(figsize = (10,8))
        plt.plot(np.arange(gTest[0]['divergeT']+1),aa[0:], 'b.')
        plt.plot(np.arange(gTest[0]['divergeT']+1), alpha[0]*np.arange(gTest[0]['divergeT']+1) + alpha[1], 'r')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    return alpha[0]

  def Dc(self,gTest,gNTests,gMaxDivergeT,gDivergence,gCSum,debug_plot=False):
    IOTA = self.IOTA
    MAX_LN_R = self.MAX_LN_R
    MIN_LN_R = self.MIN_LN_R
    N_LN_R = self.N_LN_R
    k1 = (MAX_LN_R-MIN_LN_R)/N_LN_R
    k2 = MIN_LN_R
    keepGoing = 1
    i1 = 0
    
    while keepGoing:
      for testN in range(gNTests):
        if gCSum[i1][testN]!=0:
          keepGoing = 0
          break
      i1 += keepGoing
    i1 -= 1
    if i1<0 or i1>=N_LN_R:
      i1 = 0
    keepGoing = 1
    i2 = N_LN_R-1
    while keepGoing:
      for testN in range(gNTests):
        if gCSum[i2][testN]!=0:
          keepGoing = 0
          break
      i2 -= keepGoing
    i2 += 1
    if i2<0 or i2>=N_LN_R:
      i2 = N_LN_R-1

    aa = []
    bb = []
    cc = []
    dd = []
    for i in range(i1,i2+1):
      dd.append(k1*i+k2)
      for testN in range(gNTests):
        aa.append(gCSum[i][testN])
      rs = testN
      for testN in range(rs,2*gNTests-1):
        bb.append(gCSum[i][testN])
      cc.append(gCSum[i][testN])
    
    alpha = np.polyfit(dd[0:],aa[0:],1)

    dd = np.asarray(dd)

    if debug_plot == True:
        plt.figure(figsize = (10,8))
        plt.plot(dd[0:],aa[0:], 'b.')
        plt.plot(dd[0:], alpha[0]*dd[0:] + alpha[1], 'r')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    return alpha[0]

  def compute(self,series_data,gMaxDivergeT = 600,m=5,J=7,W=100,divergeT=300,debug_lyap_plot=False,debug_Dc_plot=False):
      IOTA = self.IOTA
      MAX_LN_R = self.MAX_LN_R
      MIN_LN_R = self.MIN_LN_R
      N_LN_R = self.N_LN_R
      Test = {'filename':[],'startIndex':[],'stopIndex':[],'seriesN':series_data,'m':m,'J':J,'W':W,'divergeT':divergeT}
      gNTests = 1
      gTest = []
      gTest.append(Test)
      for i in range(gNTests):
        if gTest[i]['divergeT'] > gMaxDivergeT:
          gMaxDivergeT = 1+gTest[i]['divergeT']
      gNDivergence = np.zeros((gMaxDivergeT,1))[:,0]
      gDivergence = self.AllocateDMatrix(int(gMaxDivergeT), int(2*gNTests))
      gCSum = self.AllocateDMatrix(N_LN_R, int(2*gNTests))

      for i in range(gNTests):
        gDivergence, gCSum = self.ProcessTest(i,gTest,gMaxDivergeT,gNDivergence,gDivergence,gCSum)

      gDivergence = self.ComputeSlopes(gMaxDivergeT,gDivergence,gNTests,gCSum)

      LLE = self.lyap_r(gTest,gNTests,gMaxDivergeT,gDivergence,debug_plot=debug_lyap_plot)

      Dc_coeff = self.Dc(gTest,gNTests,gMaxDivergeT,gDivergence,gCSum,debug_plot=debug_Dc_plot)

      return LLE, Dc_coeff