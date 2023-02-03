import numpy as np
import time
from joblib import Parallel, delayed
import pickle, os
from OtherFunctions import Output

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
    
#    @staticmethod
#    def Flow(Inlet,Sim,IsoThermal): ########################### This Could be Improved "for i,j in self.Inlet.InletNodeIndexes: self.Sim[i][j].State = self.Inlet.State"
#        for Index in Inlet.InletNodeIndexes:
#            Sim[Index].State = Inlet.State
#    #        print(Index,Inlet.InletLinks)
#            Sim[Index].Gas_Particles[:] = 0# * self.Sim[Index].Gas_Particles
#            Sim[Index].Gas_Particles[:,Inlet.InletLinks] = Inlet.FlowRate
#            if IsoThermal == False:
#                Sim[Index].Energies[:] = 0# * self.Sim[Index[0]][Index[1]].Energies
#                Sim[Index].Energies[Inlet.InletLinks] = Inlet.Energy
#        return Sim
    
    @staticmethod
    def Propagation(Sim,W,H,PropagationList,Bin2Dec,ZeroState):    
        NewSim = [Cell() for i in range(W*H)]
        
        for i in range(H*W):
            NewSim[i].Type = Sim[i].Type
            NewSim[i].Solid_Particles = Sim[i].Solid_Particles
            NewSim[i].FreeVolume = Sim[i].FreeVolume
            NewSim[i].Gas_Particles[:] = 0# * self.Sim[i][j].Gas_Particles
            NewStateStr = ZeroState
            PropList = PropagationList[NewSim[i].Type]
            PropCount = len(PropList[0])####################PropList.shape[1]
            
            for k in range(PropCount):
                #print(k, PropCount, i+PropList[2][k], j+PropList[3][k])
                n_c = PropList[0][k]
                o_c = PropList[1][k]
                o_i = i+PropList[2][k]
#                    print("i={},j={},oi={},oj={},nc={},oc={}, type={}".format(i,j,o_i,o_j,n_c,o_c,NewSim[i][j].Type))
                NewSim[i].Gas_Particles[:,n_c] = Sim[o_i].Gas_Particles[:,o_c]
                NewStateStr[n_c] = Bin2Dec[Sim[o_i].State][o_c]
            if Sim[i].Energies != []:
                for k in range(PropCount):
                    NewSim[i].Energies[:,n_c] = Sim[o_i].Energies[:,o_c]
            NewStateStr = ''.join(NewStateStr)
            NewSim[i].State = Bin2Dec[NewStateStr]
        return NewSim
    
    @staticmethod
    def StaticCollision(cell, FluidCellTypes, CollisionRules, OutputCouns, randcount):
#        for j in range(self.W):
        if cell.Type in FluidCellTypes: # Internal Fluid Cell
            if cell.State in CollisionRules:
                Value = CollisionRules[cell.State]
                choice = np.random.randint(1,OutputCouns)
                #  0:Collision
                randcount = randcount+1 #Just for test
                NewState = Value[2*choice+1]
                OldIndex = Value[1]
                NewIndex = Value[2*choice]
                cell.Apply_Collision(NewState, OldIndex, NewIndex)
        return cell, randcount
                      
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
        self.randcount = np.zeros([13,Iter.Max+1], dtype = int) ######## Just for test
        self.TotalSolidWeight = TotalSolidWeight
#        self.SolidWeight = SolidWeight
#        self.FreeVolume = FreeVolume
        self.ScalingFactors = ScalingFactors
        self.ZeroState = ['0']*Neighborhood_Count
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
    
    def Flow(self):
        for Index in self.Inlet.InletNodeIndexes:
            self.Sim[Index].State = self.Inlet.State
            self.Sim[Index].Gas_Particles[:] = 0# * self.Sim[Index].Gas_Particles
            self.Sim[Index].Gas_Particles[:,self.Inlet.InletLinks] = self.Inlet.FlowRate
            if self.IsoThermal == False:
                self.Sim[Index].Energies[:] = 0# * self.Sim[Index[0]][Index[1]].Energies
                self.Sim[Index].Energies[self.Inlet.InletLinks] = self.Inlet.Energy
                
    def Collision(self):
        OutputCouns = list(self.CollisionRules.values())[0][0]
        RC = np.zeros(self.H*self.W, dtype = int)
        for i in range(self.H*self.W):
            self.Sim[i], RC[i] = self.StaticCollision(self.Sim[i], self.FluidCellTypes, self.CollisionRules, OutputCouns, 0)
        self.randcount[0,self.Iter.Current] = self.randcount[0,self.Iter.Current]+sum(RC) #Just for test
#        par_list = Parallel(n_jobs=1, prefer="threads")(delayed(self.StaticCollision)(self.Sim[i], self.FluidCellTypes, self.CollisionRules, OutputCouns, 0) for i in range(self.H*self.W))
#        self.randcount[0,self.Iter.Current] = self.randcount[0,self.Iter.Current]+sum([item[1] for item in par_list]) #Just for test
#        self.Sim = [item[0] for item in par_list]
        
#            pool = multiprocessing.Pool(4)
#            zip(*pool.map(self.Par_Col, range(0, 10 * offset, offset)))
        
    
    def Absorbtion(self, SolidComp, GS, i, k, FreeVolume):
        
        AdsCo = np.dot(SolidComp,self.AdsCoefficient[GS,:])
        n = self.Sim[i].Gas_Particles[GS,k]
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
        n = self.Sim[i].Gas_Particles[GS,0]
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
        self.Sim[i].Gas_Particles[GS,0] = self.Sim[i].Gas_Particles[GS,0] + MovingIn
        self.Sim[i].Gas_Particles[GS,k] = self.Sim[i].Gas_Particles[GS,k] - MovingIn       
#                        if MovingIn!=0:
#                            print("Adsorption: i={}, j={}, k={}, GS={}, Gas(0)={}, Gas(k)={}".format(i,j,k,GS,self.Sim[i][j].Gas_Particles[GS,0],self.Sim[i][j].Gas_Particles[GS,k]))
#                            print("            F.V.B.={}, A.C.A={}, M.In={}, F.V.A={}".format(FreeVolume,self.AdsCapacityAccupy[GS],MovingIn,FreeVolume - MovingIn*self.AdsCapacityAccupy[GS]))
#                        if self.Sim[i][j].Gas_Particles[GS,0]<0 or self.Sim[i][j].Gas_Particles[GS,k]<0:
#                            print("Adsorption: i = {}, j = {}, Gas(0) = {}, Gas(k) = {}, Move = {}".format(i,j,self.Sim[i][j].Gas_Particles[GS,0],self.Sim[i][j].Gas_Particles[GS,k],MovingIn))
        FreeVolume = FreeVolume - MovingIn*self.AdsCapacityAccupy[GS]
        self.Adsorbed[GS,self.Iter.Current] = self.Adsorbed[GS,self.Iter.Current] + MovingIn
        
        return FreeVolume
    
    
    def Difff(self, ExistingGasInLinks, k, GasOrder, FreeVolume, SolidComp, i, NeighborsList):
        if ExistingGasInLinks[k] != 0: #Adsorption And Desoprption
            for GS in GasOrder:
                FreeVolume = self.Absorbtion(SolidComp, GS, i, k, FreeVolume)
                
        i2 = i + NeighborsList[2][k]
        if i2 in self.ObstacleIndex: #Diffusion
            FreeVolume2 = self.Sim[i2].FreeVolume ############np.dot(self.Sim[i2][j2].Solid_Particles , self.AdsCapacity) - np.dot(self.Sim[i2][j2].Gas_Particles[:,0],self.AdsCapacityAccupy)
            for GS in GasOrder:
                DiffOutCo = np.dot(SolidComp,self.DiffCoefficient[GS,:])
                n = self.Sim[i].Gas_Particles[GS,0]
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
                self.Sim[i].Gas_Particles[GS,0] = self.Sim[i].Gas_Particles[GS,0] - DiffusingOut
                self.Sim[i2].Gas_Particles[GS,0] = self.Sim[i2].Gas_Particles[GS,0] + DiffusingOut            
#                        if MovingIn!=0:
#                            print("Diffusion: i={}, j={}, Gas(0,i,j)={}, i2={}, j2={}, Gas(0,i2,j2)={}".format(i,j,self.Sim[i][j].Gas_Particles[GS,0],i2,j2,self.Sim[i2][j2].Gas_Particles[GS,0]))
#                            print("           C1 : F.V.B.={}, A.C.A={}, M.In={}, F.V.A={}".format(FreeVolume,self.AdsCapacityAccupy[GS],MovingIn,FreeVolume - MovingIn*self.AdsCapacityAccupy[GS]))
#                            print("           C2 : F.V.B.={}, A.C.A={}, M.In={}, F.V.A={}".format(FreeVolume2,self.AdsCapacityAccupy[GS],MovingIn,FreeVolume2 + MovingIn*self.AdsCapacityAccupy[GS]))
#                        if self.Sim[i][j].Gas_Particles[GS,0]<0 or self.Sim[i2][j2].Gas_Particles[GS,0]<0:
#                            print("Diffusion: i = {}, j = {}, Gas(0,i,j) = {}, i = {}, j = {}, Gas(0,i2,j2) = {}, Move = {}".format(i,j,self.Sim[i][j].Gas_Particles[GS,0],i2,j2,self.Sim[i2][j2].Gas_Particles[GS,0],DiffusingOut))
                FreeVolume = FreeVolume + DiffusingOut*self.AdsCapacityAccupy[GS]
                FreeVolume2 = FreeVolume2 - DiffusingOut*self.AdsCapacityAccupy[GS]
            self.Sim[i2].FreeVolume = FreeVolume2
            
        self.Sim[i].FreeVolume = FreeVolume    
    
    def Difff_Absorbtion(self, i, NeighborOrder, GasOrder):
#        i = self.ObstacleIndex[ind]
        NeighborsList = self.PropagationList[self.Sim[i].Type]
        ExistingGasInLinks = self.Sim[i].Calculate_Link_Total_Gas_Particles()
        #            AdsorbedGasInCell = ExistingGasInLinks[0]
        FreeVolume = self.Sim[i].FreeVolume
        #            Available4Desorption = 100000 # I should add a max value for gas 
        #                # links and compute the courrunt and available volume based 
        #                # on existing gas particles and their Molar Volumes
        SolidComp = self.Sim[i].Solid_Particles/sum(self.Sim[i].Solid_Particles)
        for k in NeighborOrder:
            self.Difff(ExistingGasInLinks, k, GasOrder, FreeVolume, SolidComp, i, NeighborsList)
        
    
    def Diffusion(self):
        NeighborOrder = np.random.permutation(self.Neighborhood_Count-1)+1 # To have a 
                    #random order of diffusion and adsorption between neighbors
        GasOrder = np.random.permutation(self.Gas_Sub_Count) # To have a 
                    #random order of diffusion and adsorption between gas substances
        self.Adsorbed[:,self.Iter.Current] = self.Adsorbed[:,self.Iter.Current-1]
        
        
        #for ind in range(len(self.ObstacleIndex)):
#            print("i = {}, j = {}, Solids = {}, Gas = {}, AdsCapacity = {}, AdsCapacityAccupy = {}".format(i,j,self.Sim[i][j].Solid_Particles,self.Sim[i][j].Gas_Particles,self.AdsCapacity,self.AdsCapacityAccupy))
            #self.Difff_Absorbtion(ind, NeighborOrder, GasOrder)
        Parallel(n_jobs=1, prefer="threads")(delayed(self.Difff_Absorbtion)(ind, NeighborOrder, GasOrder) for ind in self.ObstacleIndex)#range(len(self.ObstacleIndex)))
            
        
    def Reaction(self):
        pass
            
    #0:Propagation, 1:Collision, 2:Diffusion, 3:Total
    def Simulate(self):
        while self.Iter.Current < self.Iter.Max:        
            self.Iter.Current = self.Iter.Current+1
            
#            self.Sim = self.Flow(self.Inlet,self.Sim,self.IsoThermal)
            self.Flow()
            
            tic = time.time()
            self.Sim = self.Propagation(self.Sim,self.W,self.H,self.PropagationList,self.Bin2Dec,self.ZeroState)
            self.RunTime[0,self.Iter.Current] = time.time() - tic

            tic2 = time.time()
            self.Collision()
            self.RunTime[1,self.Iter.Current] = time.time() - tic2
            
            tic2 = time.time()
            self.Diffusion()
            self.RunTime[2,self.Iter.Current] = time.time() - tic2
            
#            self.Reaction()
            
            self.RunTime[3,self.Iter.Current] = time.time() - tic
            if self.Iter.Current % self.Iter.Display == 0:
                Output(self.Iter, self.H, self.W, self.GrainSize, self.Sim, 
                       self.Gas_Sub_Count, self.Neighborhood_Count, self.LinkSpeed,
                       self.FluidCellTypes, self.ScalingFactors, self.SaveParameters, 
                       self.Adsorbed, self.RunTime, self.randcount)
            
            if self.Iter.Current % self.Iter.Save == 0:
                if os.path.isdir(self.SaveParameters.Where) == False:
                    os.mkdir(self.SaveParameters.Where)
                FullPath = self.SaveParameters.Where+'/'+self.SaveParameters.Name+str(self.Iter.Current)+'.pkl'
                pickle.dump(self, open(FullPath, "wb"))