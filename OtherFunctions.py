import pickle, os
#import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
    

def sub2ind(array_width, row, col):
    return row*array_width + col

def Set_Initial_Test_Sim(H, W, Sim):
    
    Sim[0].Type = 5#   5. Upper/Left Bundary Cell
    Sim[W-1].Type = 6#   6. Upper/Right Bundary Cell
    
    for j in range(1,W-1):
        Sim[j].Type = 3#   3. Upper/Inner Bundary Cell
        
        for i in range(2,H-1,2):
            Sim[sub2ind(W,i,j)].Type = 1#   1. Inner Fluid Cell on Even Numbered Rows
        for i in range(1,H-1,2):
            Sim[sub2ind(W,i,j)].Type = 2#   2. Inner Fluid Cell on Odd Numbered Rows
        
        Sim[sub2ind(W,H-1,j)].Type = 4#   4. Lower/Inner Bundary Cell
    
    Sim[sub2ind(W,H-1,0)].Type = 7#   7. Lower/Left Bundary Cell
    Sim[sub2ind(W,H-1,W-1)].Type = 8#   8. Lower/Right Bundary Cell
    
    for i in range(2,H-1,2):
        Sim[sub2ind(W,i,0)].Type = 10#  10. Left Entree Fluid Cell on Even Numbered Rows
        Sim[sub2ind(W,i,W-1)].Type = 12#  12. Right Inner Fluid Cell on Even Numbered Rows
    for i in range(1,H-1,2):
        Sim[sub2ind(W,i,0)].Type = 9#  9. Left Inner Fluid Cell on Odd Numbered Rows
        Sim[sub2ind(W,i,W-1)].Type = 11#  11. Right Exit Fluid Cell on Odd Numbered Rows
    
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

def OneDimPropagationTable(PropagationTable, W):
    NewTable = {}
    for k in PropagationTable:
        Value = PropagationTable[k]
        NewValue = copy.deepcopy(Value[0:3])
        for c in range(len(Value[0])):
            NewValue[2][c] = Value[2][c]*W+Value[3][c]
        NewTable[k] = NewValue
    return NewTable
    
def SquarePorousObstacle(Sim,i1,j1,Shift,Height,Width,NumColumns,BlockSize,BlockDistance,SolidCount,AdsCapacity,Mr2):#,Porosity,NumParticles,NumRows):
    
    ObstacleIndex = []
#    SolidWeight = []
#    FreeVolume = []
    TotalSolidWeight = 0
    ColShift = Shift*(int(NumColumns/2)+1)
    for b_i in range(NumColumns):
        for b_j in range(NumColumns):
            i0 = ColShift[b_j]+b_i*(BlockDistance+BlockSize)#+BlockDistance
            j0 = b_j*(BlockDistance+BlockSize)+BlockDistance
            for i in range(i0,i0+BlockSize):
                for j in range(j0,j0+BlockSize):
                    Ind = sub2ind(Width,i,j)
#                    print(i,j,Ind)
                    ObstacleIndex.append(Ind)
                    if i % 2 == 1:
                        Sim[Ind].Type = 13
                    else:
                        Sim[Ind].Type = 14
                    Sim[Ind].Solid_Particles = SolidCount
                    Sim[Ind].FreeVolume = np.dot(SolidCount , AdsCapacity)
#                    SolidWeight.append(np.dot(SolidCount , Mr2))
#                    FreeVolume.append(
                    TotalSolidWeight = TotalSolidWeight + np.dot(SolidCount , Mr2)
            
#    Col = 0
#    for q in range(j1+BlockDistance,j1+Width-BlockDistance,BlockSize+BlockDistance):
#        for k in range(i1+ColShift[Col],i1+Height-ColShift[Col+1],BlockSize+BlockDistance):
#            for i in range(k,k+BlockSize):                
#                for j in range(q,q+BlockSize):
#                    Ind = sub2ind(Width,i,j)
#                    print(i,j,Ind)
#                    ObstacleIndex.append(Ind)
#                    if i % 2 == 1:
#                        Sim[Ind].Type = 13
#                    else:
#                        Sim[Ind].Type = 14
#                    Sim[Ind].Solid_Particles = SolidCount
#                    Sim[Ind].FreeVolume = np.dot(SolidCount , AdsCapacity)
##                    SolidWeight.append(np.dot(SolidCount , Mr2))
##                    FreeVolume.append(
#                    TotalSolidWeight = TotalSolidWeight + np.dot(SolidCount , Mr2)
#        Col = Col + 1
            
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

def PlotObstacle(H,W,Sim,ADD,Name):
    Obs = np.zeros([H,W], dtype=int)
    for i in range(H):
        for j in range(W):
            Ind = sub2ind(W,i,j)
            Obs[i,j] = Sim[Ind].Type
    if os.path.isdir(ADD) == False:
        os.mkdir(ADD)
    plt.imshow(Obs)
    plt.colorbar()
    FullPath = ADD+'/Obstacle-'+Name+'.png'
    plt.savefig(FullPath, dpi=300)           
    plt.close()
    
def Output(Iter, H, W, GS, Sim, Gas_Sub_Count, Neighborhood_Count, LinkSpeed,
           FluidCellTypes, ScalingFactors, SaveParameters, Adsorbed, RunTime, randcount):
    print("Iteration = ", Iter.Current)                    
    
    GH = int(H/GS)
    GW = int(W/GS)
    TotalGasParticles = np.zeros([2,GH,GW],dtype=int)
    GasRatio = np.zeros([2*Gas_Sub_Count,GH,GW])
    TotalSolidWeight = np.zeros([GH,GW])
    ConcentrationMacro = np.zeros([2*Gas_Sub_Count,GH,GW])
    Speed = np.zeros([2,GH,GW])
    
    for i1 in range(GH):
        for j1 in range(GW):
            FluidCellCount = 0
            SolidCellCount = 0
            for i in range(i1*GS,(i1+1)*GS):
                for j in range(j1*GS,(j1+1)*GS):
                    ind = i*W + j
                    #print("i1 = %d, j1 = %d, i = %d , j = %d " %(i1,j1,i,j))
                    TempSpeed = Sim[ind].Calculate_Speed(Neighborhood_Count, LinkSpeed)
                    TempTotalGasParticles = Sim[ind].Calculate_Total_Gas_Particles()
                    Speed[0][i1][j1] = Speed[0][i1][j1]+TempSpeed[0]*TempTotalGasParticles
                    Speed[1][i1][j1] = Speed[1][i1][j1]+TempSpeed[1]*TempTotalGasParticles
                    if Sim[ind].Type in FluidCellTypes:
                        FluidCellCount = FluidCellCount + 1
                        TotalGasParticles[0,i1,j1] = TotalGasParticles[0,i1,j1] + TempTotalGasParticles
                        for k in range(Gas_Sub_Count):
                            GasRatio[2*k,i1,j1] = GasRatio[2*k,i1,j1] + Sim[ind].Calculate_Total_Sub_Gas_Particles(k)
                    else:
                        SolidCellCount = SolidCellCount + 1
                        TotalSolidWeight[i1,j1] = TotalSolidWeight[i1,j1] + np.dot(Sim[ind].Solid_Particles , ScalingFactors.Mr2)
                        Link0TotalGasParticles = Sim[ind].Calculate_Link_Total_Gas_Particles()[0]
                        TotalGasParticles[1,i1,j1] = TotalGasParticles[1,i1,j1] + Link0TotalGasParticles
                        TotalGasParticles[0,i1,j1] = TotalGasParticles[0,i1,j1] + TempTotalGasParticles - Link0TotalGasParticles                         
                        for k in range(Gas_Sub_Count):
                            GasRatio[2*k+1,i1,j1] = GasRatio[2*k+1,i1,j1] + Sim[ind].Gas_Particles[k,0]
                    

            if FluidCellCount != 0:
                if TotalGasParticles[0,i1,j1] != 0:
                    Speed[0][i1][j1] = Speed[0][i1][j1]/TotalGasParticles[0,i1,j1]
                    Speed[1][i1][j1] = Speed[1][i1][j1]/TotalGasParticles[0,i1,j1]
                    for k in range(Gas_Sub_Count):
                        ConcentrationMacro[2*k,i1,j1] = GasRatio[2*k,i1,j1]*ScalingFactors.Mr1/(FluidCellCount*ScalingFactors.MV)
                        GasRatio[2*k,i1,j1] = GasRatio[2*k,i1,j1]/TotalGasParticles[0,i1,j1]
                        
            if SolidCellCount != 0:
                if TotalGasParticles[1,i1,j1] != 0:
                    for k in range(Gas_Sub_Count):
                        GasRatio[2*k+1,i1,j1] = GasRatio[2*k+1,i1,j1]/TotalGasParticles[1,i1,j1]
                        ConcentrationMacro[2*k+1,i1,j1] = TotalGasParticles[1,i1,j1]*ScalingFactors.Mr3/TotalSolidWeight[i1,j1]
                    
    if os.path.isdir(SaveParameters.Where) == False:
        os.mkdir(SaveParameters.Where)
    
    
    if Gas_Sub_Count>1:
        cmap = "cool"
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.5))
        fig.suptitle("Total Gas Particles' Count", fontsize=16)
        images = []
        images.append(axes[0].imshow(TotalGasParticles[0,:,:], cmap=cmap))
        axes[0].set_title("Outside the Solid",fontsize=10)
        images.append(axes[1].imshow(TotalGasParticles[1,:,:], cmap=cmap))
        axes[1].set_title("Inside the Solid", fontsize=10)
        
        # Find the min and max of all colors for use in setting the color scale.
        vmin = min(image.get_array().min() for image in images)
        vmax = max(image.get_array().max() for image in images)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        for im in images:
            im.set_norm(norm)                    
        fig.colorbar(images[0], ax=axes, orientation='horizontal', fraction=.15)
        FullPath = SaveParameters.Where+'/Total Particles'+str(Iter.Current)+'.png'
        fig.savefig(FullPath, dpi=300)           
        plt.close(fig)
    
    
#            for k in range(Gas_Sub_Count):
#                fig, axes = plt.subplots(1,2)#, figsize=(5, 12))
#                fig.suptitle("Gas Ratio of the Gas Number "+str(k+1), fontsize=16)
#                images = []
#                images.append(axes[0].imshow(GasRatio[2*k,:,:]))
#                axes[0].set_title("Outside the Solid")
#                images.append(axes[1].imshow(GasRatio[2*k+1,:,:]))
#                axes[1].set_title("Inside the Solid")
#                
#                vmin = min(0,min(image.get_array().min() for image in images))
#                vmax = max(1,max(image.get_array().max() for image in images))
#                norm = colors.Normalize(vmin=vmin, vmax=vmax)
#                for im in images:
#                    im.set_norm(norm)
#                
#                fig.colorbar(images[0], ax=axes, orientation='horizontal', fraction=.15)
        
    for k in range(Gas_Sub_Count):
        fig, axes = plt.subplots(1,2, figsize=(12, 6.5))
        fig.suptitle("Gas Ratio of the Gas Number "+str(k+1), fontsize=16)
        im = axes[0].imshow(ConcentrationMacro[2*k,:,:])
#                norm = colors.Normalize(vmin=0, vmax=1)
#                im.set_norm(norm)
        axes[0].set_title("Outside the Solid", fontsize=12)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        divider = make_axes_locatable(axes[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, ax=axes[0], cax=cax)
        cbar.set_label('Adsorbate Concentration (mol/m3)', rotation=90, fontsize=10)
        axes[0].axis('off')
        cbar.formatter.set_powerlimits((-2, 3))
        cbar.update_ticks()
        im = axes[1].imshow(ConcentrationMacro[2*k+1,:,:])
#                norm = colors.Normalize(vmin=0, vmax=???)
#                im.set_norm(norm)
        axes[1].set_title("Inside the Solid", fontsize=12)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        divider = make_axes_locatable(axes[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, ax=axes[1], cax=cax)
        cbar.set_label('Adsorbed Amount (Kg/Kg)', rotation=90, fontsize=10)
        axes[1].axis('off')
        cbar.formatter.set_powerlimits((-2, 3))
        cbar.update_ticks()
        plt.tight_layout()
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)
        FullPath = SaveParameters.Where+'/Gas Ratio Gas '+str(k+1)+' '+str(Iter.Current)+'.png'
        fig.savefig(FullPath, dpi=300)
        plt.close(fig)
    
    TotalGasParticles2 = np.zeros([2,H,W],dtype=int)
    for i in range(H):
        for j in range(W):
            ind = i*W + j
            TempTotalGasParticles = Sim[ind].Calculate_Total_Gas_Particles()
            if Sim[ind].Type in FluidCellTypes:
                TotalGasParticles2[0,i,j] = TempTotalGasParticles
            else:
                Link0TotalGasParticles = Sim[ind].Calculate_Link_Total_Gas_Particles()[0]
                TotalGasParticles2[1,i,j] = Link0TotalGasParticles
                TotalGasParticles2[0,i,j] = TempTotalGasParticles - Link0TotalGasParticles
                
    cmap = "cool"
    fig, axes = plt.subplots(1,2, figsize=(12, 6.5))
    fig.suptitle("Total Gas Particles' Count", fontsize=16)
    im = axes[0].imshow(TotalGasParticles2[0,:,:])
#                norm = colors.Normalize(vmin=0, vmax=???)
#                im.set_norm(norm)
    axes[0].set_title("Outside the Solid", fontsize=12)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=axes[0], cax=cax)
    axes[0].axis('off')
    cbar.formatter.set_powerlimits((-2, 3))
    cbar.update_ticks()
    im = axes[1].imshow(TotalGasParticles2[1,:,:])
#                norm = colors.Normalize(vmin=0, vmax=???)
#                im.set_norm(norm)
    axes[1].set_title("Inside the Solid", fontsize=12)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=axes[1], cax=cax)
    axes[1].axis('off')
    cbar.formatter.set_powerlimits((-2, 3))
    cbar.update_ticks()
    plt.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)
    FullPath = SaveParameters.Where+'/Total Gas Particles\' Count '+str(Iter.Current)+'.png'
    fig.savefig(FullPath, dpi=300)
    plt.close(fig)
        
#            Y, X = np.mgrid[0:GH, 0:GW]
#            fig, ax = plt.subplots()
#            ax.streamplot(X, Y, Speed[0,:,:], Speed[1,:,:])
#            #plt.show()
#            FullPath = SaveParameters.Where+'\\Streamlines '+str(Iter.Current)+'.png'
#            fig.savefig(FullPath, format='png', dpi=300)
#            plt.close(fig)
    
    
#            for i in range(Gas_Sub_Count):
#                plt.plot(Adsorbed[i,1:Iter.Current],label = 'Gas #{}'.format(i+1))
#            plt.legend()
#            plt.xlabel('Iteration')
#            plt.ylabel('Adsorption (Particles)')
#            FullPath = SaveParameters.Where+'\\Adsorption_Micro '+str(Iter.Current)+'.png'
#            plt.savefig(FullPath, dpi=300)
#            plt.close()
    
    Time = np.array(range(1,Iter.Current))*ScalingFactors.tr*1e6
    for i in range(Gas_Sub_Count):
        plt.plot(Time,Adsorbed[i,1:Iter.Current]*ScalingFactors.Mr3[i]/TotalSolidWeight,label = 'Gas #{}'.format(i+1))
#                print(Adsorbed[i,1]*ScalingFactors.Mr3[i]/TotalSolidWeight)
#                print(Adsorbed[i,Iter.Current]*ScalingFactors.Mr3[i]/TotalSolidWeight)
    plt.legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(-2, 3))
    plt.xlabel('Time (\u03BCs)')
    plt.ylabel('Adsorption Kg/Kg')
#            FullPath = SaveParameters.Where+'/Adsorption '+str(Iter.Current)+'.png'
    FullPath = SaveParameters.Where+'/Adsorption.png'
    plt.savefig(FullPath, dpi=300)
    plt.close()
    
    
    fig, ax = plt.subplots(2, 2)
    fig.suptitle("Code Execution Time", fontsize=10)
    ax[0,0].set_title('Propagation', fontsize=6)
    ax[0,0].plot(RunTime[0,1:Iter.Current])
    ax[0,0].set_xlabel('Iteration', fontsize=6)
    ax[0,0].set_ylabel('Time (s)', fontsize=6)
    ax[0,1].set_title('Collision', fontsize=6)
    ax[0,1].plot(RunTime[1,1:Iter.Current])
    ax[0,1].set_xlabel('Iteration', fontsize=6)
    ax[0,1].set_ylabel('Time (s)', fontsize=6)
    ax[1,0].set_title('Diffusion', fontsize=6)
    ax[1,0].plot(RunTime[2,1:Iter.Current])
    ax[1,0].set_xlabel('Iteration', fontsize=6)
    ax[1,0].set_ylabel('Time (s)', fontsize=6)
    ax[1,1].set_title('Total', fontsize=6)
    ax[1,1].plot(RunTime[3,1:Iter.Current])
    ax[1,1].set_xlabel('Iteration', fontsize=6)
    ax[1,1].set_ylabel('Time (s)', fontsize=6)
    plt.tight_layout()
#            FullPath = SaveParameters.Where+'/Execution Time '+str(Iter.Current)+'.png'
    FullPath = SaveParameters.Where+'/Execution Time.png'
    fig.savefig(FullPath, format='png', dpi=300)
    plt.close(fig)
    
    
    fig, ax = plt.subplots(4, 5)
    fig.suptitle("Times of Random Number Generation", fontsize=10)
    ax[0,0].set_title('Collision', fontsize=6)
    ax[0,0].plot(randcount[0,1:Iter.Current])
    
    ax[1,0].set_title('Needed', fontsize=6)
    ax[1,0].plot(randcount[4,1:Iter.Current])
    ax[1,0].set_ylabel('Adsorption', fontsize=6)
    ax[1,1].set_title('Generated', fontsize=6)
    ax[1,1].plot(randcount[1,1:Iter.Current]+randcount[2,1:Iter.Current]+randcount[3,1:Iter.Current])
    ax[1,2].set_title('Gausian', fontsize=6)
    ax[1,2].plot(randcount[1,1:Iter.Current])
    ax[1,3].set_title('RandBank', fontsize=6)
    ax[1,3].plot(randcount[2,1:Iter.Current])
    ax[1,4].set_title('Particles', fontsize=6)
    ax[1,4].plot(randcount[3,1:Iter.Current])
    
    ax[2,0].set_title('Needed', fontsize=6)
    ax[2,0].plot(randcount[8,1:Iter.Current])
    ax[2,0].set_ylabel('Desorption', fontsize=6)
    ax[2,1].set_title('Generated', fontsize=6)
    ax[2,1].plot(randcount[5,1:Iter.Current]+randcount[6,1:Iter.Current]+randcount[7,1:Iter.Current])
    ax[2,2].set_title('Gausian', fontsize=6)
    ax[2,2].plot(randcount[5,1:Iter.Current])
    ax[2,3].set_title('RandBank', fontsize=6)
    ax[2,3].plot(randcount[6,1:Iter.Current])
    ax[2,4].set_title('Particles', fontsize=6)
    ax[2,4].plot(randcount[7,1:Iter.Current])
    
    ax[3,0].set_title('Needed', fontsize=6)
    ax[3,0].plot(randcount[12,1:Iter.Current])
    ax[3,0].set_ylabel('Diffusion', fontsize=6)
    ax[3,1].set_title('Generated', fontsize=6)
    ax[3,1].plot(randcount[9,1:Iter.Current]+randcount[10,1:Iter.Current]+randcount[11,1:Iter.Current])
    ax[3,2].set_title('Gausian', fontsize=6)
    ax[3,2].plot(randcount[9,1:Iter.Current])
    ax[3,3].set_title('RandBank', fontsize=6)
    ax[3,3].plot(randcount[10,1:Iter.Current])
    ax[3,4].set_title('Particles', fontsize=6)
    ax[3,4].plot(randcount[11,1:Iter.Current])
    plt.tight_layout()
#            FullPath = SaveParameters.Where+'/Random Generation '+str(Iter.Current)+'.png'
    FullPath = SaveParameters.Where+'/Random Generation.png'
    fig.savefig(FullPath, format='png', dpi=300)
    plt.close(fig)