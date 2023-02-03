from inspect import currentframe, getframeinfo
import sys
#----Should be removed at the end------
import os
from Classes import Cell, Iteration, Save, Inlet_Set, Lattice_Set, Scaling
from Dictionaries import PropagationTable
from OtherFunctions import Set_Initial_Test_Sim, SquarePorousObstacle
from OtherFunctions import SingleRoundObstacle, SingleSquareObstacle
from OtherFunctions import QSGS, OneDimPropagationTable# ,MakeCollisionRules
from OtherFunctions import PlotObstacle
import pickle
import numpy as np


def Initialize(StartChoice=[], SimParams=[], ObstacleChoice=[], ObstacleParams=[]):
    
    if StartChoice == []:
        print("How do you want to start the simulation?")
        print("- L: Load from previous simulation file")
        print("- N: New Simulation")
        print("- C: Changed minde (Cancel)")
        StartChoice = input("Enter L / N / C : ")
    
    if ((StartChoice == 'N') or (StartChoice == 'n')):
        Neighborhood_Count=7
        LinkSpeed = [np.array([0,0]),np.array([1,0]),
                     np.array([0.5,-np.sqrt(3)/2]),
                     np.array([-0.5,-np.sqrt(3)/2]),
                     np.array([-1,0]),np.array([-0.5,np.sqrt(3)/2]),
                     np.array([0.5,np.sqrt(3)/2])]
        Cs = 340.4
        cs = np.sqrt(3/7)
        ur = Cs/cs            
            
        if SimParams == []:
            MH = int(input('Enter Macroscopic Height : ')) # Macroscopic Value of the H
            MW = int(input('Enter Macroscopic Width : '))   # Macroscopic Value of the W
            H = int(input('Enter Microscopic Height : ')) # Microscopic Value of the H
            Gas_Sub_Count = int(input('Enter Number of Gas Substances : '))
            Flow = np.zeros(Gas_Sub_Count, dtype='int') # Number of Gas Particles in Enterance Cells
            for i in range(Gas_Sub_Count):
                Flow[i] = int(input('Enter Number of Gas {} Particles in Enterance : '.format(i+1)))
            Solid_Sub_Count = int(input('Enter Number of Solid Substances : '))
            SolidCount = np.zeros(Solid_Sub_Count, dtype='int') # Number of Gas Particles in Enterance Cells
            for i in range(Solid_Sub_Count):
                SolidCount[i] = int(input('Enter Number of Solid {} Particles in Each Cell : '.format(i+1)))
            IsoThermal = True
            GrainSize = int(input('Enter Grain Size : '))
            
            MT = float(input('Enter Macroscopic Simulation Time : '))  # Macroscopic Value of the Final Simulation Time
            SaveIter = int(input('Enter Save iteration number : '))
            DispIter = int(input('Enter Display iteration number : '))
            
            GMM = float(input('Enter Gas Molar Mass (Kg/mol) : ')) # Gas Molar Mass (Kg/mol)
            GD = float(input('Enter Gas Adsorbate Density (Kg/m^3) : ')) # Gas Adsorbate Density (Kg/m^3)
            MZ = float(input('Enter Thickness of the Cell (m) : ')) # Macroscopic Thickness of the Cell (m)
            
            M_AdsCoefficient = np.array([[1.47e7]]) # Adsorption rate constant (m3 mol-1 s-1) [Solid_Sub_Count,Gas_Sub_Count]
            M_DesCoefficient = np.array([[2.94e6]]) # Desorption rate constant (s-1) [Solid_Sub_Count,Gas_Sub_Count]
            M_DiffCoefficient = np.array([[3.92e-8]]) # Solid diffusion coefficient (m^2/s) Solid_Sub_Count,Gas_Sub_Count]
            
            M_AdsCapacity = np.array([0.32]) #Saturation adsorption Weight (Kg/Kg) [Solid_Sub_Count,Gas_Sub_Count]
            SD = np.array([670]) # Solid Density (Kg/m3)
            GD = np.array([0.58]) # Gas Density (Kg/m3)
            
            Energy = []
            
            ADD = input('Enter folder name to save data : ')
            Name = input('Enter file name to save data : ')        
            
        else:
            MH = SimParams[0] # Macroscopic Value of the H
            H = SimParams[1] ####################int(input('Enter Height : '))
            MW = SimParams[2]   # Macroscopic Value of the W
            
            Gas_Sub_Count = SimParams[3]
            Flow = SimParams[4]
            Solid_Sub_Count = SimParams[5]
            SolidCount = SimParams[6]
            IsoThermal = SimParams[7]
            GrainSize = SimParams[8]
            
            MT = SimParams[9]  # Macroscopic Value of the Final Simulation Time
            SaveIter = SimParams[10] # Save iteration number
            DispIter = SimParams[11] # Display iteration number
            
            GMM = SimParams[12] # Gas Molar Mass (Kg/mol)
            GD = SimParams[13] # Gas Adsorbate Density (Kg/m^3)
            MZ = SimParams[14] # Macroscopic Thickness of the Cell (m)
            
            M_AdsCoefficient = SimParams[15] # Adsorption rate constant (m3 mol-1 s-1) [Solid_Sub_Count,Gas_Sub_Count]
            M_DesCoefficient = SimParams[16] # Desorption rate constant (s-1) [Solid_Sub_Count,Gas_Sub_Count]
            M_DiffCoefficient = SimParams[17] # Solid diffusion coefficient (m^2/s) Solid_Sub_Count,Gas_Sub_Count]
            
            M_AdsCapacity = SimParams[18] #Saturation adsorption Weight (Kg/Kg) [Solid_Sub_Count,Gas_Sub_Count]
            SD = SimParams[19] # Solid Density (Kg/m3)
            GD = SimParams[20] # Gas Density (Kg/m3)
            
            Energy = SimParams[21]
            
            ADD = SimParams[22] # Folder name to save data
            Name = SimParams[23] # File name to save data        
        
        Lr = MH/H # Scaling factor for length
        W = int(MW/Lr)#########int(input('Enter Width : '))
        tr = Lr/ur # Scaling factor for time  
        MaxIter = int(MT/tr) ####################int(input('Enter final iteration number : '))
        Iter = Iteration(0, SaveIter, DispIter, MaxIter)
        
        MGPC = sum(Flow)*7 # Maximum Gas Particles in each cell
        MSPC = sum(SolidCount) # Maximum Solid Particles in each cell
        MA = 1.5*np.sqrt(3)*Lr**2  # Area of Each Cell (m^2)
        MV = MZ*MA # Volume of Each Cell (m^3)
        MGC = GD*MV # Maximum Gas in each cell (Kg)
        PG1M = MGPC*GMM/MGC # Number of Gas Particles for 1 mol of Gas
        Mr1 = 1/PG1M
        
        AdsCoefficient = M_AdsCoefficient*Mr1*tr/Lr**3
        DesCoefficient = M_DesCoefficient*tr
        DiffCoefficient = M_DiffCoefficient*tr/Lr**2
                          
        MSWC =  SD*MV # (Kg) Maximum Solid Weight in each cell
        MGWC =  GD*MV # (Kg) Maximum Gas Weight in each cell
        PS1K = MSPC/MSWC # Number of Solid Particles for 1 Kg of Solid
        PG1K = MGPC/MGWC # Number of Gas Particles for 1 Kg of Gas
        Mr2 = 1/PS1K
        Mr3 = 1/PG1K
        AdsCapacity = int(M_AdsCapacity*Mr2/Mr3)
        AdsCapacityAccupy = np.array([1]) #[Gas_Sub_Count]
        
        ScalingFactors = Scaling(Lr,tr,Mr1,Mr2,Mr3,MV)
        
       
        InletNodeIndexes = list(range(2,H-1,2))
        InletNodeIndexes = [W*x for x in InletNodeIndexes]
        Inlet = Inlet_Set(InletNodeIndexes, [1,2,6], Flow, Energy,70)
        
        CurrentDirectory = os.getcwd()
        ADD = CurrentDirectory+'/'+ADD
        SaveParameters = Save(ADD, Name)
        
        Sim = [Cell(Neighborhood_Count, Gas_Sub_Count, Solid_Sub_Count, IsoThermal) for i in range(W*H)]
        Sim = Set_Initial_Test_Sim(H, W, Sim)

#        Sim, ObstacleIndex = Set_Initial_Test_Sim(H, W, Sim)
        if ObstacleChoice == []:
            print("How do you want to design the obstacle?")
            print("- 0: No Obstacle")
            print("- 1: Load from previous simulation file")
            print("- 2: SquarePorousObstacle")
            print("- 3: SingleRoundObstacle")
            print("- 4: SingleSquareObstacle")
            print("- 5: QSGS")
            ObstacleChoice = int(input("Enter 0~5 : "))
            
        if (ObstacleChoice == 0):
            ObstacleIndex = []
            TotalSolidWeight = 0
        elif (ObstacleChoice == 1):
            cf = currentframe()
            filename = getframeinfo(cf).filename
            print("\n\nObstacle Load should be completed at\n\t",filename,"\n\t Line Number : ",cf.f_lineno)
            
        elif (ObstacleChoice == 2):        
            print("Designing Porous Media Similar to LBM...")
            
            if ObstacleParams == []:
                i1 = int(input("Enter Obstacle\'s 1st Row Number : "))
                j1 = int(input("Enter Obstacle\'s 1st Column Number : "))
                Height = int(input("Enter Obstacle\'s Height : "))
                Width = int(input("Enter Obstacle\'s Width : "))
                NumColumns = int(input("Enter Number of Columns/Blocks : "))
                BlockSize = int(input("Enter Block Size : "))
                BlockDistance = int(input("Enter Obstacle\'s 1st Row number : "))
                Shift = [int(1.5*BlockDistance),int(0.5*BlockDistance)]#########BlockDistance#int(BlockSize/2)
            else:
                i1 = ObstacleParams[0]
                j1 = ObstacleParams[1]
                Height = ObstacleParams[2]
                Width = ObstacleParams[3]
                NumColumns = ObstacleParams[4]
                BlockSize = ObstacleParams[5]
                BlockDistance = ObstacleParams[6]
                Shift = ObstacleParams[7]
                
            Sim, ObstacleIndex, TotalSolidWeight = SquarePorousObstacle(Sim,i1,j1,Shift,Height,Width,NumColumns,BlockSize,BlockDistance,SolidCount,AdsCapacity,Mr2)
            PlotObstacle(Height,Width,Sim,ADD,ObstacleParams[8])
        elif (ObstacleChoice == 3):
            ic = int(H/2)
            jc = int(W/5)
            Diameter = int(2*H/5)
            Sim, ObstacleIndex = SingleRoundObstacle(Sim,ic,jc,Diameter,SolidCount,AdsCapacity)
        
        elif (ObstacleChoice == 4):
            i1 = int(3*H/10)
            j1 = int(W/5)
            Height = int(4*H/10)
            Width = Height
            Sim, ObstacleIndex = SingleSquareObstacle(Sim,i1,j1,Height,Width,SolidCount,AdsCapacity)
        
        elif (ObstacleChoice == 5):
            ic = int(H/2)
            jc = int(W/5)
            Diameter = int(4.5*H/5)
            Porosity = 0.6
            Cd = 0.02
            Di = 0.5
            Sim, ObstacleIndex = QSGS(Sim,ic,jc,Diameter,Porosity,Cd,Di,SolidCount,AdsCapacity)
        
        
#        PropagationList = pickle.load(open('PropagationList.p',"rb"))
        PropagationList = OneDimPropagationTable(PropagationTable,W)
#        NeighborCount = Sim[0][0].Gas_Particles.shape[1]
#        CollisionRules = MakeCollisionRules(NeighborCount)
        CollisionRules = pickle.load(open('CollisionRules.pkl',"rb"))
        Bin2Dec = pickle.load(open('Bin2DecDic.pkl',"rb"))
        
        FluidCellTypes = [1,2,9,10,11,12] #####################################
        SolidCellTypes = [13,14] #####################################
        
        RandBankLength = int(1e6)
        
        if os.path.isdir(ADD) == False:
            os.mkdir(ADD)
        fp = open(ADD+"/Parameters.txt","w") 
        
        fp.write("Height : Macroscopit = {} , Microscopit = {} , Scaling Factor (Lr) = {}".format(MH,H,Lr))
        fp.write("\nTime : Macroscopit = {} , Microscopit = {} , Scaling Factor (tr) = {}".format(MT,MaxIter,tr))
        fp.write("\nScaling Factors : Mr1 = {} , Mr2 = {} , Mr3 = {}".format(Mr1,Mr2,Mr3))
        fp.write("\nSaturation Adsorption capacity : Macroscopit = {} , Microscopit = {}".format(M_AdsCapacity,AdsCapacity))
        fp.write("\nAdsorption rate constant : Macroscopit = {} , Microscopit = {}".format(M_AdsCoefficient,AdsCoefficient))
        fp.write("\nDesorption rate constant : Macroscopit = {} , Microscopit = {}".format(M_DesCoefficient,DesCoefficient))
        fp.write("\nSolid diffusion coefficient : Macroscopit = {} , Microscopit = {}".format(M_DiffCoefficient,DiffCoefficient))
        
        fp.write("\nPropagation List:\n")
        for kk in PropagationList:
            fp.write("\n{} : {}".format(kk,PropagationList[kk]))
        fp.write("\nCollision Rules:\n")
        for kk in CollisionRules:
            fp.write("\n{} : {}".format(kk,CollisionRules[kk]))
        fp.write("\nBin2Dec:\n")
        for kk in Bin2Dec:
            fp.write("\n{} : {}".format(kk,Bin2Dec[kk]))
        
        Lattice = Lattice_Set(H, W, Sim, Iter, Inlet, PropagationList,
                              CollisionRules, Bin2Dec, FluidCellTypes, 
                              SolidCellTypes, True, SaveParameters, LinkSpeed, 
                              ScalingFactors, GrainSize, ObstacleIndex, Neighborhood_Count,
                              Gas_Sub_Count, Solid_Sub_Count, AdsCoefficient,
                              AdsCapacity, AdsCapacityAccupy, DesCoefficient,
                              DiffCoefficient,RandBankLength,[],TotalSolidWeight)#,SolidWeight,FreeVolume)
        
    elif ((StartChoice == 'L') or (StartChoice == 'l')):
        #----Should be removed at the end------
        cf = currentframe()
        filename = getframeinfo(cf).filename
        print("\n\nLoad Simulation should be completed at\n\t",filename,"\n\t Line Number : ",cf.f_lineno)
        #----Should be removed at the end------
        ADD = 'out'#############################input('Enter folder name to load data : ')
        ADD = CurrentDirectory+'/'+ADD
        Name = 'Data200.pkl'##################input('Enter file name to load data : ')
        FullPath = ADD+'/'+Name
        Lattice = pickle.load(open(FullPath,"rb"))
        fp = open(ADD+"/Parameters.txt","a+") 
        fp.write("\n\nContinue from file : "+Name)
#        Lattice.Iter.Max = 5000
        Lattice.Iter.Display = 1000
        Lattice.Iter.Save = 1000
        if len(Lattice.Adsorbed[0])<=Lattice.Iter.Max:
            Lattice.Adsorbed = np.append(Lattice.Adsorbed, np.zeros([Lattice.Gas_Sub_Count,Lattice.Iter.Max-len(Lattice.Adsorbed[0])+1]), axis=1)
    elif ((StartChoice == 'C') or (StartChoice == 'c')):
        print("Exiting simulation...\n\n")
        #----Should be removed at the end------
        cf = currentframe()
        filename = getframeinfo(cf).filename
        print("\n\nBetter implementing of exit (without error) should be replaced at \n\t",filename,"\n\t Line Number : ",cf.f_lineno)
        #----Should be removed at the end------
        sys.exit()
    else:
        print("Wrong choice, try again...\n\n")
        Lattice = Initialize();
        
    fp.close() 
    return Lattice