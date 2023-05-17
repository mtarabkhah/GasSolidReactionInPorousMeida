import pickle
#import matplotlib.pyplot as plt
import numpy as np

def cls():
    print ("\n"*80)
    
def TempPrint(txt):
    from inspect import currentframe, getframeinfo
    cf = currentframe()
    filename = getframeinfo(cf).filename
    print("{} \n\t {}\n\t Line Number : {}".format(txt,filename, cf.f_lineno))
    
    
def MyBinary(i,NeighborCount):
    bin_i = np.zeros(NeighborCount, dtype=int)
    for k in range(NeighborCount):
        bin_i[NeighborCount-k-1] = i%2
        i = int(i/2)
    return bin_i

def MyBinary2(n,NeighborCount):
    bin_i = MyBinary(n,NeighborCount)
    IdexList = []
    for k in range(len(bin_i)):
        if bin_i[k] != 0:
            IdexList = IdexList + [NeighborCount-k-1]
    return IdexList
    
def MyBinaryString(n,NeighborCount):
    bin_i = MyBinary(n,NeighborCount)
    BinStr = ''
    for i in range(NeighborCount-1,-1,-1):
        BinStr = BinStr + str(bin_i[i])
    return BinStr

def MakeBin2DecDic():
    Bin2Dec = {}
    NeighborCount = 7
    for i in range(2**NeighborCount):
        BinStr = MyBinaryString(i,NeighborCount)
        Bin2Dec[i] = BinStr
        Bin2Dec[BinStr] = i
        
    pickle.dump(Bin2Dec, open('Bin2DecDic.pkl', "wb"))
    print(Bin2Dec)
    
#MakeBin2DecDic()
    
def MakeCollisionRules2(NeighborCount,CollisionTable):
#    from Dictionaries import CollisionTable
    
    Pos = []
    for k in CollisionTable:
        Value = CollisionTable[k]
        OldIndex = [Value[1]]
        NewValue = [Value[0]] + OldIndex
        for choice in range(Value[0]):
#            OldIndex = np.concatenate((Value[1],Value[2*(choice+1)])) 
#            NewIndex = np.concatenate((Value[2*(choice+1)],Value[1]))
            
            NewIndex = Value[2*(choice+1)]
            NewValue = NewValue + [NewIndex] + [Value[2*(choice+1)+1]]
        CollisionTable[k] = NewValue
        if not(Value[0] in Pos):
            Pos.append(Value[0])
    MaxPos = max(Pos)
    
    CollisionRules = {}
    for i in CollisionTable:
        Value = CollisionTable[i]
        CollisionRules[i] = [MaxPos]+ [Value[1]] + int(MaxPos/Value[0])*Value[2:]
#    for i in range(2**NeighborCount):
#        if i in CollisionTable:
#            Value = CollisionTable[i]
#            CollisionRules[i] = [MaxPos]+int(MaxPos/Value[0])*Value[1:]
#        else:
##            Binary_i = MyBinary(i,NeighborCount)
##            Value = [MaxPos]+MaxPos*[Binary_i,Binary_i,i]
#            Value = [MaxPos]+MaxPos*[[0],[0],i]
#            CollisionRules[i] = Value
    
    pickle.dump(CollisionRules, open('CollisionRules.pkl', "wb"))
    return CollisionRules

    
def MyDec(List):
    Dec = 0
    for k in List:
        Dec = Dec+2**k
        
    return Dec
    
def MakeCollisionRules1(NeighborCount):
    from Dictionaries import CollisionTable1
    
    CollisionTable = {}
    for k in CollisionTable1:
        Value = CollisionTable1[k]
        InitialState = MyDec(Value[0])
        Count = len(Value)-1
        NewValue = [Count, np.array(Value[0])]
        for i in range(1,Count+1):
            NewValue = NewValue + [np.array(Value[i]),MyDec(Value[i])]
        CollisionTable[InitialState] = NewValue
        
    for k in range(1,2**NeighborCount):
        if not(k in CollisionTable):
            InitialState = k
            Count = 1
            OldIndex = MyBinary2(k,NeighborCount)
            NewValue = [Count, OldIndex]
            NewIndex = []
            for i in OldIndex:
                if i==0:
                    NewIndex.append(i)
                elif i<=int((NeighborCount)/2):
                    NewIndex.append(i+int((NeighborCount)/2))
                else:
                    NewIndex.append(i-int((NeighborCount)/2))
            NewValue = NewValue + [NewIndex,MyDec(NewIndex)]
            CollisionTable[InitialState] = NewValue
        
#    for k in CollisionTable:
#        print(k,':',CollisionTable[k])
        
    CollisionRules = MakeCollisionRules2(NeighborCount,CollisionTable)
    pickle.dump(CollisionRules, open('CollisionRules.pkl', "wb"))
    return CollisionRules

#
#CollisionRules = MakeCollisionRules1(7)
#for k in CollisionRules:
#    print(k,':',CollisionRules[k])
    
    
def Set_Initial_Test_Sim(H, W, Sim):
    
    Sim[0][0].Type = 5#   5. Upper/Left Bundary Cell
    Sim[0][W-1].Type = 6#   6. Upper/Right Bundary Cell
    
    for j in range(1,W-1):
        Sim[0][j].Type = 3#   3. Upper/Inner Bundary Cell
        
        for i in range(2,H-1,2):
            Sim[i][j].Type = 1#   1. Inner Fluid Cell on Even Numbered Rows
        for i in range(1,H-1,2):
            Sim[i][j].Type = 2#   2. Inner Fluid Cell on Odd Numbered Rows
        
        Sim[H-1][j].Type = 4#   4. Lower/Inner Bundary Cell
    
    Sim[H-1][0].Type = 7#   7. Lower/Left Bundary Cell
    Sim[H-1][W-1].Type = 8#   8. Lower/Right Bundary Cell
    
    for i in range(2,H-1,2):
        Sim[i][0].Type = 10#  10. Left Entree Fluid Cell on Even Numbered Rows
        Sim[i][W-1].Type = 12#  12. Right Inner Fluid Cell on Even Numbered Rows
    for i in range(1,H-1,2):
        Sim[i][0].Type = 9#  9. Left Inner Fluid Cell on Odd Numbered Rows
        Sim[i][W-1].Type = 11#  11. Right Exit Fluid Cell on Odd Numbered Rows
    
    return Sim
#    ObstacleIndex = []
#    for i in range(int(H/3),int(2*H/3)):
#        for j in range(int(W/5),int(W/5)+1):#int(W/5),int(2*W/5)):
#            ObstacleIndex.append([i,j])
#            if i % 2 == 1:
#                Sim[i][j].Type = 13
#            else:
#                Sim[i][j].Type = 14
#            Sim[i][j].Solid_Particles = 1
#    return Sim, ObstacleIndex

def SquarePorousObstacle(Sim,i1,j1,Shift,Height,Width,NumColumns,BlockSize,BlockDistance,SolidCount,AdsCapacity,Mr2):#,Porosity,NumParticles,NumRows):
    
    ObstacleIndex = []
#    SolidWeight = []
#    FreeVolume = []
    TotalSolidWeight = 0
    ColShift = Shift*(int(NumColumns/2)+1)
    Col = 0
    for q in range(j1+BlockDistance,j1+Width-BlockDistance,BlockSize+BlockDistance):
        for k in range(i1+ColShift[Col],i1+Height-ColShift[Col+1],BlockSize+BlockDistance):
            for i in range(k,k+BlockSize):
                for j in range(q,q+BlockSize):
                    ObstacleIndex.append([i,j])
                    if i % 2 == 1:
                        Sim[i][j].Type = 13
                    else:
                        Sim[i][j].Type = 14
                    Sim[i][j].Solid_Particles = SolidCount
                    Sim[i][j].FreeVolume = np.dot(SolidCount , AdsCapacity)
#                    SolidWeight.append(np.dot(SolidCount , Mr2))
#                    FreeVolume.append(
                    TotalSolidWeight = TotalSolidWeight + np.dot(SolidCount , Mr2)
        Col = Col + 1
            
    return Sim, ObstacleIndex, TotalSolidWeight#, SolidWeight, FreeVolume

def SingleSquareObstacle(Sim,i1,j1,Height,Width,SolidCount,AdsCapacity):
    
    ObstacleIndex = []
    for i in range(i1,i1+Height):
        for j in range(j1,j1+Width):
            ObstacleIndex.append([i,j])
            if i % 2 == 1:
                Sim[i][j].Type = 13
            else:
                Sim[i][j].Type = 14
            Sim[i][j].Solid_Particles = SolidCount
            Sim[i][j].FreeVolume = np.dot(SolidCount , AdsCapacity)
            
    return Sim, ObstacleIndex

def SingleRoundObstacle(Sim,ic,jc,Diameter,SolidCount,AdsCapacity):
    
    Radius = int(Diameter/2)
    ObstacleIndex = []
    for i in range(ic-Radius-1,ic+Radius+1):
        for j in range(jc-Radius-1,jc+Radius+1):
            if np.sqrt((i-ic)**2+(j-jc)**2) <= Radius:
                ObstacleIndex.append([i,j])
                if i % 2 == 1:
                    Sim[i][j].Type = 13
                else:
                    Sim[i][j].Type = 14
                Sim[i][j].Solid_Particles = SolidCount
                Sim[i][j].FreeVolume = np.dot(SolidCount , AdsCapacity)
            
    return Sim, ObstacleIndex

def QSGS(Sim,ic,jc,Diameter,Porosity,Cd,Di,SolidCount,AdsCapacity):
    
    Radius = int(Diameter/2)
    TotalObstacleIndex = []
    for i in range(ic-Radius-1,ic+Radius+1):
        for j in range(jc-Radius-1,jc+Radius+1):
            if np.sqrt((i-ic)**2+(j-jc)**2) <= Radius:
                TotalObstacleIndex.append([i,j])    
    TotalCellCount = len(TotalObstacleIndex)
    
    r = np.random.rand(TotalCellCount)
    CoreInd = np.argwhere(r<Cd)
    ObstacleIndex = []
    for k in range (TotalCellCount-1,-1,-1):
        if k in CoreInd:
            i,j = TotalObstacleIndex[k]
            ObstacleIndex.append([i,j])     
    SolidCellCount = len(ObstacleIndex)
    while (SolidCellCount/TotalCellCount) < (1-Porosity):
        for k in range(SolidCellCount):
            i1,j1 = ObstacleIndex[k]
            for i in range(i1-1,i1+2):
                for j in range(j1-1,j1+2):
                    if ([i,j] in TotalObstacleIndex) and not([i,j] in ObstacleIndex):
                        ObstacleIndex.append([i,j])
        SolidCellCount = len(ObstacleIndex)
    
    while (SolidCellCount/TotalCellCount) > (1-Porosity):
        del ObstacleIndex[-1]
        SolidCellCount = SolidCellCount - 1
        
    for i,j in ObstacleIndex:
        if i % 2 == 1:
            Sim[i][j].Type = 13
        else:
            Sim[i][j].Type = 14
        Sim[i][j].Solid_Particles = SolidCount
        Sim[i][j].FreeVolume = np.dot(SolidCount , AdsCapacity)
        
    print("final porosity = %s"%(1-(SolidCellCount/TotalCellCount)))
    return Sim, ObstacleIndex