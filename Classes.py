import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pickle, os
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
        

class Cell:
    def __init__(self, Neighborhood_Count=7, Gas_Sub_Count=1, Solid_Sub_Count=1, IsoThermal=True, Type=1, State=0, FreeVolume = 1):
        self.Type = Type
        self.Gas_Particles = np.zeros([Gas_Sub_Count,Neighborhood_Count], dtype=int)
        self.Solid_Particles = np.zeros(Solid_Sub_Count, dtype=int)
        if IsoThermal:            
            self.Energies = []
        else:
            self.Energies = np.zeros(Neighborhood_Count, dtype=int)
        self.State = State
        self.FreeVolume = FreeVolume
#        self.Volume = self.CalculateVolume(Molar_Volume)
        
    def Calculate_Volume(self,Molar_Volume):
        return sum(self.Solid_Particles*Molar_Volume)
    
    def Define_State(self,Neighborhood_Count):
#        Dec = 0        
#        print(np.nonzero(np.count_nonzero(self.Gas_Particles, axis=0)))
#        for k in np.nonzero(np.count_nonzero(self.Gas_Particles, axis=0)):
#            Dec = Dec+2**k
#        print(Dec)
#        self.State = int(Dec)
        
        Dec = 0
        for k in range(Neighborhood_Count):#len(self.Gas_Particles[0])):
            if any(self.Gas_Particles[:,k] != 0):
                Dec = Dec+2**k
        self.State = Dec
        
    def Apply_Collision(self, State, OldIndex, NewIndex):
        self.State = State
        Temp = self.Gas_Particles[:,OldIndex]
        self.Gas_Particles[:,OldIndex] = 0
        self.Gas_Particles[:,NewIndex] = Temp
        if self.Energies != []:
            Temp = self.Energies[OldIndex]
            self.Energies[OldIndex] = 0
            self.Energies[NewIndex] = Temp
    
    def Calculate_Total_Gas_Particles(self):
        return sum(sum(self.Gas_Particles))
    
    def Calculate_Link_Total_Gas_Particles(self):
        return sum(self.Gas_Particles)
    
    def Calculate_Total_Sub_Gas_Particles(self,Substance):
        return sum(self.Gas_Particles[Substance,:])
    
    def Calculate_GasRatio(self, Substance):
        return sum(self.Gas_Particles[Substance,:])/sum(sum(self.Gas_Particles))
    
    def Calculate_Speed(self, LinkCount, LinkSpeed):
        Speed = np.zeros(2)
        ParticleCount = self.Calculate_Link_Total_Gas_Particles()
        for i in range(LinkCount):
            Speed = Speed + ParticleCount[i]*LinkSpeed[i][:]        
        return Speed
    
    # Cell types:
    #   1. Inner Fluid Cell on Odd Numbered Rows
    #   2. Inner Fluid Cell on Even Numbered Rows
    #   3. Upper/Inner Bundary Cell
    #   4. Lower/Inner Bundary Cell
    #   5. Upper/Left Bundary Cell
    #   6. Upper/Right Bundary Cell
    #   7. Lower/Left Bundary Cell
    #   8. Lower/Right Bundary Cell
    #   9. Left Inner Fluid Cell on Odd Numbered Rows
    #  10. Left Entree Fluid Cell on Even Numbered Rows
    #  11. Right Exit Fluid Cell on Odd Numbered Rows
    #  12. Right Inner Fluid Cell on Even Numbered Rows
    #  13. Obstacle Cell on Odd Numbered Rows
    #  14. Obstacle Cell on Even Numbered Rows
    
class Iteration:
    def __init__(self, Current=0, Save=100, Display=200, Max=1000):
        self.Current = Current
        self.Save = Save
        self.Display = Display
        self.Max = Max
        
class Save:
    def __init__(self,Where='Results', Name='Data'):
        self.Where = Where
        self.Name = Name
        
class Inlet_Set:
    def __init__(self,InletNodeIndexes, InletLinks=[1,2,6], FlowRate=1, Energy=[],State=70):
        self.InletNodeIndexes = InletNodeIndexes
        self.InletLinks = InletLinks
        self.FlowRate = FlowRate
        self.Energy = Energy
        self.State = State

class Scaling:
    def __init__(self,Lr,tr,Mr1,Mr2,Mr3,MV):
        self.Lr = Lr
        self.tr = tr
        self.Mr1 = Mr1
        self.Mr2 = Mr2
        self.Mr3 = Mr3
        self.MV = MV
        
class Lattice_Set:
    def __init__(self, H, W, Sim, Iter, Inlet, PropagationList, CollisionRules,
                 Bin2Dec, FluidCellTypes, SolidCellTypes, IsoThermal,
                 SaveParameters, LinkSpeed, ScalingFactors, GrainSize = 8, ObstacleIndex =[],
                 Neighborhood_Count=7, Gas_Sub_Count=1, Solid_Sub_Count=1,
                 AdsCoefficient = 0, AdsCapacity = 0, AdsCapacityAccupy = 0,
                 DesCoefficient = 0, DiffCoefficient = 0,RandBankLength = 1e6, Adsorbed = [],TotalSolidWeight = 0):#,SolidWeight = 0,FreeVolume = 0):
        self.H = H
        self.W = W
        self.Sim = Sim
        self.Iter = Iter
        self.SaveParameters = SaveParameters
        self.Inlet = Inlet
        self.PropagationList = PropagationList
        self.CollisionRules = CollisionRules
        self.FluidCellTypes = FluidCellTypes
        self.SolidCellTypes = SolidCellTypes
        self.IsoThermal = IsoThermal
        self.Neighborhood_Count = Neighborhood_Count
        self.Gas_Sub_Count = Gas_Sub_Count
        self.Solid_Sub_Count = Solid_Sub_Count
        self.ObstacleIndex = ObstacleIndex
        self.GrainSize = GrainSize
        self.LinkSpeed = LinkSpeed
        self.Bin2Dec = Bin2Dec
        self.AdsCoefficient = AdsCoefficient
        self.AdsCapacity = AdsCapacity
        self.AdsCapacityAccupy = AdsCapacityAccupy
        self.DesCoefficient = DesCoefficient
        self.DiffCoefficient = DiffCoefficient
        self.RandBank = np.random.rand(RandBankLength)
        self.randcount = np.zeros([13,Iter.Max+1]) ######## Just for test
        self.TotalSolidWeight = TotalSolidWeight
#        self.SolidWeight = SolidWeight
#        self.FreeVolume = FreeVolume
        self.ScalingFactors = ScalingFactors
        #  0:Collision
        #  1:Adsorption_Gausian
        #  2:Adsorption_RandBank
        #  3:Adsorption_Particles
        #  4:Adsorption_Needed
        #  5:Desorption_Gausian
        #  6:Desorption_RandBank
        #  7:Desorption_Particles
        #  8:Desorption_Needed
        #  9:Diffusion_Gausian
        # 10:Diffusion_RandBank
        # 11:Diffusion_Particles
        # 12:Diffusion_Needed
        self.RunTime = np.zeros([4,Iter.Max+1]) ######## Just for test
        #0:Propagation, 1:Collision, 2:Diffusion, 3:Total
        if Adsorbed == []:
            self.Adsorbed = np.zeros([Gas_Sub_Count,Iter.Max+1], dtype=int)
        elif len(Adsorbed[0])<=Iter.Max:
            self.Adsorbed = np.append(Adsorbed, np.zeros([Gas_Sub_Count,Iter.Max-len(Adsorbed[0])+1], dtype=int), axis=1)
        
    def Flow(self): ########################### This Could be Improved "for i,j in self.Inlet.InletNodeIndexes: self.Sim[i][j].State = self.Inlet.State"
        for i in range(len(self.Inlet.InletNodeIndexes)):
            Index = self.Inlet.InletNodeIndexes[i]
            self.Sim[Index[0]][Index[1]].State = self.Inlet.State
    #        print(Index,Inlet.InletLinks)
            self.Sim[Index[0]][Index[1]].Gas_Particles = 0 * self.Sim[Index[0]][Index[1]].Gas_Particles
            self.Sim[Index[0]][Index[1]].Gas_Particles[:,self.Inlet.InletLinks] = self.Inlet.FlowRate
            if self.IsoThermal == False:
                self.Sim[Index[0]][Index[1]].Energies = 0 * self.Sim[Index[0]][Index[1]].Energies
                self.Sim[Index[0]][Index[1]].Energies[self.Inlet.InletLinks] = self.Inlet.Energy
                
#        for i,j in self.Inlet.InletNodeIndexes:
#            self.Sim[i][j].State = self.Inlet.State
#    #        print(Index,Inlet.InletLinks)
#            self.Sim[i][j].Gas_Particles = 0 * self.Sim[i][j].Gas_Particles
#            self.Sim[i][j].Gas_Particles[:,self.Inlet.InletLinks] = self.Inlet.FlowRate
#            if self.IsoThermal == False:
#                self.Sim[i][j].Energies = 0 * self.Sim[i][j].Energies
#                self.Sim[i][j].Energies[self.Inlet.InletLinks] = self.Inlet.Energy
        
    def Propagation(self):    
        NewSim = [[Cell() for i in range(self.W)] for j in range(self.H)]
        
        for j in range(self.W):
            for i in range(self.H):        
                NewSim[i][j].Type = self.Sim[i][j].Type
                NewSim[i][j].Solid_Particles = self.Sim[i][j].Solid_Particles
                NewSim[i][j].FreeVolume = self.Sim[i][j].FreeVolume
                NewSim[i][j].Gas_Particles = 0 * self.Sim[i][j].Gas_Particles
                NewStateStr = ['0']*self.Neighborhood_Count
                PropList = self.PropagationList[NewSim[i][j].Type]
                PropCount = len(PropList[0])####################PropList.shape[1]
                
                for k in range(PropCount):
                    #print(k, PropCount, i+PropList[2][k], j+PropList[3][k])
                    n_c = PropList[0][k]
                    o_i = i+PropList[2][k]
                    o_j = j+PropList[3][k]
                    o_c = PropList[1][k]
#                    print("i={},j={},oi={},oj={},nc={},oc={}, type={}".format(i,j,o_i,o_j,n_c,o_c,NewSim[i][j].Type))
                    NewSim[i][j].Gas_Particles[:,n_c] = self.Sim[o_i][o_j].Gas_Particles[:,o_c]
                    NewStateStr[n_c] = self.Bin2Dec[self.Sim[o_i][o_j].State][o_c]
                if self.Sim[i][j].Energies != []:
                    for k in range(PropCount):
                        NewSim[i][j].Energies[:,n_c] = self.Sim[o_i][o_j].Energies[:,o_c]
                NewStateStr = ''.join(NewStateStr)
                NewSim[i][j].State = self.Bin2Dec[NewStateStr]
        self.Sim = NewSim
        
    def Collision(self):
        OutputCouns = list(self.CollisionRules.values())[0][0]
        for i in range(self.H): #----Should be confirmed at the end------ (Just for test)
            for j in range(self.W):
                if self.Sim[i][j].Type in self.FluidCellTypes: # Internal Fluid Cell
                    if self.Sim[i][j].State in self.CollisionRules:
                        Value = self.CollisionRules[self.Sim[i][j].State]
                        choice = np.random.randint(1,OutputCouns)
                        #  0:Collision
                        self.randcount[0,self.Iter.Current] = self.randcount[0,self.Iter.Current]+1 #Just for test
                        NewState = Value[2*choice+1]
                        OldIndex = Value[1]
                        NewIndex = Value[2*choice]
                        self.Sim[i][j].Apply_Collision(NewState, OldIndex, NewIndex)
    
    
    def Diffusion(self):
        NeighborOrder = np.random.permutation(self.Neighborhood_Count-1)+1 # To have a 
                    #random order of diffusion and adsorption between neighbors
        GasOrder = np.random.permutation(self.Gas_Sub_Count) # To have a 
                    #random order of diffusion and adsorption between gas substances
        self.Adsorbed[:,self.Iter.Current] = self.Adsorbed[:,self.Iter.Current-1]
        for i,j in self.ObstacleIndex:
#            print("i = {}, j = {}, Solids = {}, Gas = {}, AdsCapacity = {}, AdsCapacityAccupy = {}".format(i,j,self.Sim[i][j].Solid_Particles,self.Sim[i][j].Gas_Particles,self.AdsCapacity,self.AdsCapacityAccupy))
            NeighborsList = self.PropagationList[self.Sim[i][j].Type]
            ExistingGasInLinks = self.Sim[i][j].Calculate_Link_Total_Gas_Particles()
#            AdsorbedGasInCell = ExistingGasInLinks[0]
            FreeVolume = self.Sim[i][j].FreeVolume
#            Available4Desorption = 100000 # I should add a max value for gas 
#                # links and compute the courrunt and available volume based 
#                # on existing gas particles and their Molar Volumes
            SolidComp = self.Sim[i][j].Solid_Particles/sum(self.Sim[i][j].Solid_Particles)
            for k in NeighborOrder:
                if ExistingGasInLinks[k] != 0: #Adsorption And Desoprption
                    for GS in GasOrder:
                        AdsCo = np.dot(SolidComp,self.AdsCoefficient[GS,:])
                        n = self.Sim[i][j].Gas_Particles[GS,k]
                        MaxAds = min(n,int(FreeVolume/self.AdsCapacityAccupy[GS]))
                        #  1:Adsorption_Gausian
                        #  2:Adsorption_RandBank
                        #  3:Adsorption_Particles
                        #  4:Adsorption_Needed
                        self.randcount[4,self.Iter.Current] = self.randcount[4,self.Iter.Current]+n #Just for test
                        Mean = n*AdsCo
                        Variance = Mean * (1-AdsCo)
                        if n<=len(self.RandBank):
                            start = np.random.randint(0,len(self.RandBank)-n)
                            Adsorbing = min(MaxAds, len(np.argwhere(self.RandBank[start:start+n]<AdsCo)))
                            self.randcount[2,self.Iter.Current] = self.randcount[2,self.Iter.Current]+1 #Just for test
                        elif Mean >= 9 and Variance>=0:
                            Adsorbing = max(0,min(np.random.normal(Mean, Variance),n))
                            self.randcount[1,self.Iter.Current] = self.randcount[1,self.Iter.Current]+1 #Just for test
                        else:
                            Adsorbing = min(MaxAds, len(np.argwhere(np.random.rand(n)<AdsCo)))
                            self.randcount[3,self.Iter.Current] = self.randcount[3,self.Iter.Current]+n #Just for test
#                        Desorbing = min(Available4Desorption, len(np.argwhere(np.random.rand(self.Sim[i][j].Gas_Particles[GS,0])<DesCo)))
                        DesCo = np.dot(SolidComp,self.DesCoefficient[GS,:])
                        n = self.Sim[i][j].Gas_Particles[GS,0]
                        #  5:Desorption_Gausian
                        #  6:Desorption_RandBank
                        #  7:Desorption_Particles
                        #  8:Desorption_Needed
                        self.randcount[8,self.Iter.Current] = self.randcount[8,self.Iter.Current]+n #Just for test
                        Mean = n*DesCo
                        Variance = Mean * (1-DesCo)
                        if n<=len(self.RandBank):
                            start = np.random.randint(0,len(self.RandBank)-n)
                            Desorbing = len(np.argwhere(self.RandBank[start:start+n]<DesCo))
                            self.randcount[6,self.Iter.Current] = self.randcount[6,self.Iter.Current]+1 #Just for test
                        elif Mean >= 9 and Variance>=0:
                            Desorbing = max(0,min(np.random.normal(Mean, Variance),n))
                            self.randcount[5,self.Iter.Current] = self.randcount[5,self.Iter.Current]+1 #Just for test
                        else:
                            Desorbing = len(np.argwhere(np.random.rand(n)<DesCo))
                            self.randcount[7,self.Iter.Current] = self.randcount[7,self.Iter.Current]+n #Just for test
                        MovingIn = Adsorbing - Desorbing
                        self.Sim[i][j].Gas_Particles[GS,0] = self.Sim[i][j].Gas_Particles[GS,0] + MovingIn
                        self.Sim[i][j].Gas_Particles[GS,k] = self.Sim[i][j].Gas_Particles[GS,k] - MovingIn       
#                        if MovingIn!=0:
#                            print("Adsorption: i={}, j={}, k={}, GS={}, Gas(0)={}, Gas(k)={}".format(i,j,k,GS,self.Sim[i][j].Gas_Particles[GS,0],self.Sim[i][j].Gas_Particles[GS,k]))
#                            print("            F.V.B.={}, A.C.A={}, M.In={}, F.V.A={}".format(FreeVolume,self.AdsCapacityAccupy[GS],MovingIn,FreeVolume - MovingIn*self.AdsCapacityAccupy[GS]))
#                        if self.Sim[i][j].Gas_Particles[GS,0]<0 or self.Sim[i][j].Gas_Particles[GS,k]<0:
#                            print("Adsorption: i = {}, j = {}, Gas(0) = {}, Gas(k) = {}, Move = {}".format(i,j,self.Sim[i][j].Gas_Particles[GS,0],self.Sim[i][j].Gas_Particles[GS,k],MovingIn))
                        FreeVolume = FreeVolume - MovingIn*self.AdsCapacityAccupy[GS]
                        self.Adsorbed[GS,self.Iter.Current] = self.Adsorbed[GS,self.Iter.Current] + MovingIn
                    
                i2 = i + NeighborsList[2][k]
                j2 = j + NeighborsList[3][k]
                if [i2,j2] in self.ObstacleIndex: #Diffusion
                    FreeVolume2 = self.Sim[i2][j2].FreeVolume ############np.dot(self.Sim[i2][j2].Solid_Particles , self.AdsCapacity) - np.dot(self.Sim[i2][j2].Gas_Particles[:,0],self.AdsCapacityAccupy)
                    for GS in GasOrder:
                        DiffOutCo = np.dot(SolidComp,self.DiffCoefficient[GS,:])
                        n = self.Sim[i][j].Gas_Particles[GS,0]
                        MaxDiffOut = min(n,int(FreeVolume2/self.AdsCapacityAccupy[GS]))
                        #  9:Diffusion_Gausian
                        # 10:Diffusion_RandBank
                        # 11:Diffusion_Particles
                        # 12:Diffusion_Needed
                        self.randcount[12,self.Iter.Current] = self.randcount[12,self.Iter.Current]+n #Just for test
                        Mean = n*DiffOutCo
                        Variance = Mean * (1-DiffOutCo)
                        if n<=len(self.RandBank):
                            start = np.random.randint(0,len(self.RandBank)-n)
                            DiffusingOut = min(MaxDiffOut, len(np.argwhere(self.RandBank[start:start+n]<DiffOutCo)))
                            self.randcount[10,self.Iter.Current] = self.randcount[10,self.Iter.Current]+1 #Just for test
                        elif Mean >= 9 and Variance>=0:
                            DiffusingOut = max(0,min(np.random.normal(Mean, Variance),n))
                            self.randcount[9,self.Iter.Current] = self.randcount[9,self.Iter.Current]+1 #Just for test
                        else:
                            DiffusingOut = min(MaxDiffOut, len(np.argwhere(np.random.rand(n)<DiffOutCo)))
                            self.randcount[11,self.Iter.Current] = self.randcount[11,self.Iter.Current]+n #Just for test
                        self.Sim[i][j].Gas_Particles[GS,0] = self.Sim[i][j].Gas_Particles[GS,0] - DiffusingOut
                        self.Sim[i2][j2].Gas_Particles[GS,0] = self.Sim[i2][j2].Gas_Particles[GS,0] + DiffusingOut            
#                        if MovingIn!=0:
#                            print("Diffusion: i={}, j={}, Gas(0,i,j)={}, i2={}, j2={}, Gas(0,i2,j2)={}".format(i,j,self.Sim[i][j].Gas_Particles[GS,0],i2,j2,self.Sim[i2][j2].Gas_Particles[GS,0]))
#                            print("           C1 : F.V.B.={}, A.C.A={}, M.In={}, F.V.A={}".format(FreeVolume,self.AdsCapacityAccupy[GS],MovingIn,FreeVolume - MovingIn*self.AdsCapacityAccupy[GS]))
#                            print("           C2 : F.V.B.={}, A.C.A={}, M.In={}, F.V.A={}".format(FreeVolume2,self.AdsCapacityAccupy[GS],MovingIn,FreeVolume2 + MovingIn*self.AdsCapacityAccupy[GS]))
#                        if self.Sim[i][j].Gas_Particles[GS,0]<0 or self.Sim[i2][j2].Gas_Particles[GS,0]<0:
#                            print("Diffusion: i = {}, j = {}, Gas(0,i,j) = {}, i = {}, j = {}, Gas(0,i2,j2) = {}, Move = {}".format(i,j,self.Sim[i][j].Gas_Particles[GS,0],i2,j2,self.Sim[i2][j2].Gas_Particles[GS,0],DiffusingOut))
                        FreeVolume = FreeVolume + DiffusingOut*self.AdsCapacityAccupy[GS]
                        FreeVolume2 = FreeVolume2 - DiffusingOut*self.AdsCapacityAccupy[GS]
                    self.Sim[i2][j2].FreeVolume = FreeVolume2
                    
                self.Sim[i][j].FreeVolume = FreeVolume
            
        
    def Reaction(self):
        pass
        
    def Output(self, mode):
        if mode == 1 or mode == 3:
            print("Iteration = ", self.Iter.Current)                    
        
        if mode == 2:
            GS = self.GrainSize
            GH = int(self.H/GS)
            GW = int(self.W/GS)
            TotalGasParticles = np.zeros([2,GH,GW],dtype=int)
            GasRatio = np.zeros([2*self.Gas_Sub_Count,GH,GW])
            TotalSolidWeight = np.zeros([GH,GW])
            ConcentrationMacro = np.zeros([2*self.Gas_Sub_Count,GH,GW])
            Speed = np.zeros([2,GH,GW])
            
            for i1 in range(GH):
                for j1 in range(GW):
                    FluidCellCount = 0
                    SolidCellCount = 0
                    for i in range(i1*GS,(i1+1)*GS):
                        for j in range(j1*GS,(j1+1)*GS):
                            #print("i1 = %d, j1 = %d, i = %d , j = %d " %(i1,j1,i,j))
                            TempSpeed = self.Sim[i][j].Calculate_Speed(self.Neighborhood_Count, self.LinkSpeed)
                            TempTotalGasParticles = self.Sim[i][j].Calculate_Total_Gas_Particles()
                            Speed[0][i1][j1] = Speed[0][i1][j1]+TempSpeed[0]*TempTotalGasParticles
                            Speed[1][i1][j1] = Speed[1][i1][j1]+TempSpeed[1]*TempTotalGasParticles
                            if self.Sim[i][j].Type in self.FluidCellTypes:
                                FluidCellCount = FluidCellCount + 1
                                TotalGasParticles[0,i1,j1] = TotalGasParticles[0,i1,j1] + TempTotalGasParticles
                                for k in range(self.Gas_Sub_Count):
                                    GasRatio[2*k,i1,j1] = GasRatio[2*k,i1,j1] + self.Sim[i][j].Calculate_Total_Sub_Gas_Particles(k)
                            else:
                                SolidCellCount = SolidCellCount + 1
                                TotalSolidWeight[i1,j1] = TotalSolidWeight[i1,j1] + np.dot(self.Sim[i][j].Solid_Particles , self.ScalingFactors.Mr2)
                                Link0TotalGasParticles = self.Sim[i][j].Calculate_Link_Total_Gas_Particles()[0]
                                TotalGasParticles[1,i1,j1] = TotalGasParticles[1,i1,j1] + Link0TotalGasParticles
                                TotalGasParticles[0,i1,j1] = TotalGasParticles[0,i1,j1] + TempTotalGasParticles - Link0TotalGasParticles                         
                                for k in range(self.Gas_Sub_Count):
                                    GasRatio[2*k+1,i1,j1] = GasRatio[2*k+1,i1,j1] + self.Sim[i][j].Gas_Particles[k,0]
                            
    
                    if FluidCellCount != 0:
                        if TotalGasParticles[0,i1,j1] != 0:
                            Speed[0][i1][j1] = Speed[0][i1][j1]/TotalGasParticles[0,i1,j1]
                            Speed[1][i1][j1] = Speed[1][i1][j1]/TotalGasParticles[0,i1,j1]
                            for k in range(self.Gas_Sub_Count):
                                ConcentrationMacro[2*k,i1,j1] = GasRatio[2*k,i1,j1]*self.ScalingFactors.Mr1/(FluidCellCount*self.ScalingFactors.MV)
                                GasRatio[2*k,i1,j1] = GasRatio[2*k,i1,j1]/TotalGasParticles[0,i1,j1]
                                
                    if SolidCellCount != 0:
                        if TotalGasParticles[1,i1,j1] != 0:
                            for k in range(self.Gas_Sub_Count):
                                GasRatio[2*k+1,i1,j1] = GasRatio[2*k+1,i1,j1]/TotalGasParticles[1,i1,j1]
                                ConcentrationMacro[2*k+1,i1,j1] = TotalGasParticles[1,i1,j1]*self.ScalingFactors.Mr3/TotalSolidWeight[i1,j1]
                            
            if os.path.isdir(self.SaveParameters.Where) == False:
                os.mkdir(self.SaveParameters.Where)
            
            
            if self.Gas_Sub_Count>1:
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
    
                plt.tight_layout()
    #            plt.show()
                FullPath = self.SaveParameters.Where+'/Total Particles'+str(self.Iter.Current)+'.png'
                fig.savefig(FullPath, dpi=300)           
                plt.close(fig)
            
            
#            for k in range(self.Gas_Sub_Count):
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
                
            for k in range(self.Gas_Sub_Count):
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
                cbar.formatter.set_powerlimits((0, 0))
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
                cbar.formatter.set_powerlimits((0, 0))
                cbar.update_ticks()
                plt.tight_layout()
#                plt.show()
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)
                FullPath = self.SaveParameters.Where+'/Gas Ratio Gas '+str(k+1)+' '+str(self.Iter.Current)+'.png'
                fig.savefig(FullPath, dpi=300)
                plt.close(fig)
                
                            
            TotalGasParticles2 = np.zeros([2,self.H,self.W],dtype=int)
            for i in range(self.H):
                for j in range(self.W):
                    TempTotalGasParticles = self.Sim[i][j].Calculate_Total_Gas_Particles()
                    if self.Sim[i][j].Type in self.FluidCellTypes:
                        TotalGasParticles2[0,i,j] = TempTotalGasParticles
                    else:
                        Link0TotalGasParticles = self.Sim[i][j].Calculate_Link_Total_Gas_Particles()[0]
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
            cbar.formatter.set_powerlimits((0, 0))
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
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()
            plt.tight_layout()
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)
            FullPath = self.SaveParameters.Where+'/Total Gas Particles\' Count '+str(self.Iter.Current)+'.png'
            fig.savefig(FullPath, dpi=300)
            plt.close(fig)
                
#            Y, X = np.mgrid[0:GH, 0:GW]
#            fig, ax = plt.subplots()
#            ax.streamplot(X, Y, Speed[0,:,:], Speed[1,:,:])
#            #plt.show()
#            FullPath = self.SaveParameters.Where+'/Streamlines '+str(self.Iter.Current)+'.png'
#            fig.savefig(FullPath, format='png', dpi=300)
#            plt.close(fig)
            
            
#            for i in range(self.Gas_Sub_Count):
#                plt.plot(self.Adsorbed[i,1:self.Iter.Current],label = 'Gas #{}'.format(i+1))
#            plt.legend()
#            plt.xlabel('Iteration')
#            plt.ylabel('Adsorption (Particles)')
#            FullPath = self.SaveParameters.Where+'/Adsorption_Micro '+str(self.Iter.Current)+'.png'
#            plt.savefig(FullPath, dpi=300)
#            plt.close()
            
            Time = np.array(range(1,self.Iter.Current))*self.ScalingFactors.tr*1e6
            for i in range(self.Gas_Sub_Count):
                plt.plot(Time,self.Adsorbed[i,1:self.Iter.Current]*self.ScalingFactors.Mr3[i]/self.TotalSolidWeight,label = 'Gas #{}'.format(i+1))
#                print(self.Adsorbed[i,1]*self.ScalingFactors.Mr3[i]/self.TotalSolidWeight)
#                print(self.Adsorbed[i,self.Iter.Current]*self.ScalingFactors.Mr3[i]/self.TotalSolidWeight)
            plt.legend()
            plt.xlabel('Time (\u03BCs)')
            plt.ylabel('Adsorption Kg/Kg')
            FullPath = self.SaveParameters.Where+'/Adsorption '+str(self.Iter.Current)+'.png'
            plt.savefig(FullPath, dpi=300)
            plt.close()
            
            
            fig, ax = plt.subplots(2, 2)
            fig.suptitle("Code Execution Time", fontsize=10)
            ax[0,0].set_title('Propagation', fontsize=6)
            ax[0,0].plot(self.RunTime[0,1:self.Iter.Current])
            ax[0,0].set_xlabel('Iteration', fontsize=6)
            ax[0,0].set_ylabel('Time (s)', fontsize=6)
            ax[0,1].set_title('Collision', fontsize=6)
            ax[0,1].plot(self.RunTime[1,1:self.Iter.Current])
            ax[0,1].set_xlabel('Iteration', fontsize=6)
            ax[0,1].set_ylabel('Time (s)', fontsize=6)
            ax[1,0].set_title('Diffusion', fontsize=6)
            ax[1,0].plot(self.RunTime[2,1:self.Iter.Current])
            ax[1,0].set_xlabel('Iteration', fontsize=6)
            ax[1,0].set_ylabel('Time (s)', fontsize=6)
            ax[1,1].set_title('Total', fontsize=6)
            ax[1,1].plot(self.RunTime[3,1:self.Iter.Current])
            ax[1,1].set_xlabel('Iteration', fontsize=6)
            ax[1,1].set_ylabel('Time (s)', fontsize=6)
            plt.tight_layout()
            FullPath = self.SaveParameters.Where+'/Execution Time '+str(self.Iter.Current)+'.png'
            fig.savefig(FullPath, format='png', dpi=300)
            plt.close(fig)
            
            
            fig, ax = plt.subplots(4, 5)
            fig.suptitle("Times of Random Number Generation", fontsize=10)
            ax[0,0].set_title('Collision', fontsize=6)
            ax[0,0].plot(self.randcount[0,1:self.Iter.Current])
            
            ax[1,0].set_title('Needed', fontsize=6)
            ax[1,0].plot(self.randcount[4,1:self.Iter.Current])
            ax[1,0].set_ylabel('Adsorption', fontsize=6)
            ax[1,1].set_title('Generated', fontsize=6)
            ax[1,1].plot(self.randcount[1,1:self.Iter.Current]+self.randcount[2,1:self.Iter.Current]+self.randcount[3,1:self.Iter.Current])
            ax[1,2].set_title('Gausian', fontsize=6)
            ax[1,2].plot(self.randcount[1,1:self.Iter.Current])
            ax[1,3].set_title('RandBank', fontsize=6)
            ax[1,3].plot(self.randcount[2,1:self.Iter.Current])
            ax[1,4].set_title('Particles', fontsize=6)
            ax[1,4].plot(self.randcount[3,1:self.Iter.Current])
            
            ax[2,0].set_title('Needed', fontsize=6)
            ax[2,0].plot(self.randcount[8,1:self.Iter.Current])
            ax[2,0].set_ylabel('Desorption', fontsize=6)
            ax[2,1].set_title('Generated', fontsize=6)
            ax[2,1].plot(self.randcount[5,1:self.Iter.Current]+self.randcount[6,1:self.Iter.Current]+self.randcount[7,1:self.Iter.Current])
            ax[2,2].set_title('Gausian', fontsize=6)
            ax[2,2].plot(self.randcount[5,1:self.Iter.Current])
            ax[2,3].set_title('RandBank', fontsize=6)
            ax[2,3].plot(self.randcount[6,1:self.Iter.Current])
            ax[2,4].set_title('Particles', fontsize=6)
            ax[2,4].plot(self.randcount[7,1:self.Iter.Current])
            
            ax[3,0].set_title('Needed', fontsize=6)
            ax[3,0].plot(self.randcount[12,1:self.Iter.Current])
            ax[3,0].set_ylabel('Diffusion', fontsize=6)
            ax[3,1].set_title('Generated', fontsize=6)
            ax[3,1].plot(self.randcount[9,1:self.Iter.Current]+self.randcount[10,1:self.Iter.Current]+self.randcount[11,1:self.Iter.Current])
            ax[3,2].set_title('Gausian', fontsize=6)
            ax[3,2].plot(self.randcount[9,1:self.Iter.Current])
            ax[3,3].set_title('RandBank', fontsize=6)
            ax[3,3].plot(self.randcount[10,1:self.Iter.Current])
            ax[3,4].set_title('Particles', fontsize=6)
            ax[3,4].plot(self.randcount[11,1:self.Iter.Current])
            plt.tight_layout()
            FullPath = self.SaveParameters.Where+'/Random Generation '+str(self.Iter.Current)+'.png'
            fig.savefig(FullPath, format='png', dpi=300)
            plt.close(fig)
                        
        if mode == 3:
            if os.path.isdir(self.SaveParameters.Where) == False:
                os.mkdir(self.SaveParameters.Where)
            FullPath = self.SaveParameters.Where+'/'+self.SaveParameters.Name+str(self.Iter.Current)+'.pkl'
            pickle.dump(self, open(FullPath, "wb"))
            
    #0:Propagation, 1:Collision, 2:Diffusion, 3:Total
    def Simulate(self):
        while self.Iter.Current < self.Iter.Max:        
            self.Iter.Current = self.Iter.Current+1
            
            self.Flow()
            tic = time.time()
            self.Propagation()
            toc = time.time()
            self.RunTime[0,self.Iter.Current] = toc - tic
#            print('propagation : ',toc - tic)
            tic2 = time.time()
            self.Collision()
            toc = time.time()
            self.RunTime[1,self.Iter.Current] = toc - tic2
#            print('collision : ',toc - tic2)
            
            tic2 = time.time()
            self.Diffusion()
            toc = time.time()
            self.RunTime[2,self.Iter.Current] = toc - tic2
            
#            self.Reaction()
            
            if self.Iter.Current % self.Iter.Display == 0:
                self.Output(2)
            
            if self.Iter.Current % self.Iter.Save == 0:
                self.Output(3)
            toc = time.time()
#            print('Iteration : ',toc - tic)
            toc = time.time()
            self.RunTime[3,self.Iter.Current] = toc - tic