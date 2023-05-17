def Initialize():
    #----Should be removed at the end------
    from inspect import currentframe, getframeinfo
    import sys
    #----Should be removed at the end------
    import os
    from Classes import Cell, Iteration, Save, Inlet_Set, Lattice_Set, Scaling
    from Dictionaries import PropagationTable
    from OtherFunctions import Set_Initial_Test_Sim, SquarePorousObstacle, SingleRoundObstacle, SingleSquareObstacle, QSGS# ,MakeCollisionRules
    import pickle
    import numpy as np


    CurrentDirectory = os.getcwd()
    
    print("Initializing the Simulation...")
    
    print("How do you want to start the simulation?")
    print("- L: Load from previous simulation file")
    print("- N: New Simulation")
    print("- C: Changed minde (Cancel)")
    choice = input("Enter L / N / C : ") ####################'l' #
    
    if (choice.lower() == 'n'):
        
        print("Setting Desired Input Values...")
        H = int(input('Enter Height : ')) #201 ####################
        MH = 4e-5 # Macroscopic Value of the H
        MW = MH   # Macroscopic Value of the W
        Lr = MH/H # Scaling factor for length
        W = int(input('Enter Width : ')) #int(MW/Lr) #401 ####################
        Neighborhood_Count=7
        LinkSpeed = [np.array([0,0]),np.array([1,0]),
                     np.array([0.5,-np.sqrt(3)/2]),
                     np.array([-0.5,-np.sqrt(3)/2]),
                     np.array([-1,0]),np.array([-0.5,np.sqrt(3)/2]),
                     np.array([0.5,np.sqrt(3)/2])]
        Gas_Sub_Count=1########2############
        Flow = np.array([3])############[[20],[10]]############
        Solid_Sub_Count=1########2############
        SolidCount = np.array([10])##########np.array([20,5])######
        IsoThermal=True
        GrainSize = 4 ##################
        
        Cs = 340.4
        cs = np.sqrt(3/7)
        ur = Cs/cs
        tr = Lr/ur # Scaling factor for time
        MT = 1.2e-5  # Macroscopic Value of the Final Simulation Time
        MaxIter = int(input('Enter final iteration number : ')) #int(MT/tr) ####################
        SaveIter = int(input('Enter Save iteration number : ')) #200 ####################
        DispIter = int(input('Enter Display iteration number : ')) #1000 ####################
        Iter = Iteration(0, SaveIter, DispIter, MaxIter)
        
        MGPC = sum(Flow)*7 # Maximum Gas Particles in each cell
        GMM = 0.018 # Gas Molar Mass (Kg/mol)
        GD = 0.58 # Gas Adsorbate Density (Kg/m^3)
        MA = 1.5*np.sqrt(3)*Lr**2  # Area of Each Cell (m^2)
        MZ = 1.6e-7 # Macroscopic Thickness of the Cell (m)
        MV = MZ*MA # Volume of Each Cell (m^3)
        MGC = GD*MV # Maximum Gas in each cell (Kg)
        PG1M = MGPC*GMM/MGC
        Mr1 = 1/PG1M
        
        M_AdsCoefficient = np.array([[1.47e7]]) # Adsorption rate constant (m3 mol-1 s-1) [Solid_Sub_Count,Gas_Sub_Count]
        M_DesCoefficient = np.array([[2.94e6]]) # Desorption rate constant (s-1) [Solid_Sub_Count,Gas_Sub_Count]
        M_DiffCoefficient = np.array([[3.92e-8]]) # Solid diffusion coefficient (m^2/s) Solid_Sub_Count,Gas_Sub_Count]
        
        AdsCoefficient = M_AdsCoefficient*Mr1*tr/Lr**3
        DesCoefficient = M_DesCoefficient*tr
        DiffCoefficient = M_DiffCoefficient*tr/Lr**2
        
        M_AdsCapacity = np.array([0.32]) #Saturation adsorption Weight (Kg/Kg) [Solid_Sub_Count,Gas_Sub_Count]
        MSPC = sum(SolidCount) # Maximum Solid Particles in each cell
        SD = np.array([670]) # (Kg/m3) Solid Density
        GD = np.array([0.58]) # (Kg/m3) Gas Density
        MSWC =  SD*MV # (Kg) Maximum Solid Weight in each cell
        MGWC =  GD*MV # (Kg) Maximum Gas Weight in each cell
        PS1K = MSPC/MSWC # Number of Solid Particles for 1 Kg of Solid
        PG1K = MGPC/MGWC # Number of Gas Particles for 1 Kg of Gas
        Mr2 = 1/PS1K
        Mr3 = 1/PG1K
        AdsCapacity = int(M_AdsCapacity*Mr2/Mr3)
        AdsCapacityAccupy = np.array([1]) #[Gas_Sub_Count]
        
        ScalingFactors = Scaling(Lr,tr,Mr1,Mr2,Mr3,MV)
        
        Energy = []
        InletNodeIndexes = np.zeros((len(list(range(2,H-1,2))),2), dtype=int)
        InletNodeIndexes[:,0] = list(range(2,H-1,2))
        Inlet = Inlet_Set(InletNodeIndexes, [1,2,6], Flow, Energy,70)
        
        ADD = input('Enter folder name to save data : ') # 'out' ####################
        ADD = CurrentDirectory+'/'+ADD
        Name = input('Enter file name to save data : ') #'Data' ####################
        SaveParameters = Save(ADD, Name)
        
        Sim = [[Cell(Neighborhood_Count, Gas_Sub_Count, Solid_Sub_Count, IsoThermal) for i in range(W)] for j in range(H)]
        Sim = Set_Initial_Test_Sim(H, W, Sim)
        
        

#        Sim, ObstacleIndex = Set_Initial_Test_Sim(H, W, Sim)
        
        print("How do you want to design the obstacle?")
        print("- 1: Load from previous simulation file")
        print("- 2: SquarePorousObstacle")
        print("- 3: SingleRoundObstacle")
        print("- 4: SingleSquareObstacle")
        print("- 5: QSGS")
        choice = input("Enter 1~5 : ") #'2' ####################
        
        if (choice == '1'):
            cf = currentframe()
            filename = getframeinfo(cf).filename
            print("\n\nObstacle Load should be completed at\n\t",filename,"\n\t Line Number : ",cf.f_lineno)
            
        elif (choice == '2'):      
            print("Designing Porous Media Similar to LBM...")
            i1 = 1#############int(2*H/10)
            j1 = 1#############int(W/5)
            Height = H############int(6*H/10)
            Width = W
            NumColumns = 8
            BlockSize = int((H-1)/10)###########int(3*H/20)-1
            BlockDistance = int(1.2*BlockSize)
            Shift = [int(1.8*BlockSize),int(0.6*BlockSize)]#########BlockDistance#int(BlockSize/2)
            Sim, ObstacleIndex, TotalSolidWeight = SquarePorousObstacle(Sim,i1,j1,Shift,Height,Width,NumColumns,BlockSize,BlockDistance,SolidCount,AdsCapacity,Mr2)
        
        elif (choice == '3'):
            ic = int(H/2)
            jc = int(W/5)
            Diameter = int(2*H/5)
            Sim, ObstacleIndex = SingleRoundObstacle(Sim,ic,jc,Diameter,SolidCount,AdsCapacity)
        
        elif (choice == '4'):
            i1 = int(3*H/10)
            j1 = int(W/5)
            Height = int(4*H/10)
            Width = Height
            Sim, ObstacleIndex = SingleSquareObstacle(Sim,i1,j1,Height,Width,SolidCount,AdsCapacity)
        
        elif (choice == '5'):
            ic = int(H/2)
            jc = int(W/5)
            Diameter = int(4.5*H/5)
            Porosity = 0.6
            Cd = 0.02
            Di = 0.5
            Sim, ObstacleIndex = QSGS(Sim,ic,jc,Diameter,Porosity,Cd,Di,SolidCount,AdsCapacity)
        
        
#        PropagationList = pickle.load(open('PropagationList.p',"rb"))
        PropagationList = PropagationTable
#        NeighborCount = Sim[0][0].Gas_Particles.shape[1]
#        CollisionRules = MakeCollisionRules(NeighborCount)
        CollisionRules = pickle.load(open('CollisionRules.pkl',"rb"))
        Bin2Dec = pickle.load(open('Bin2DecDic.pkl',"rb"))
        
        FluidCellTypes = [1,2,9,10,11,12] #####################################
        SolidCellTypes = [13,14] #####################################
        
        Smallest = 1e-6##################min(min())
        RandBankLength = int(9/Smallest)
        
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
        
        Lattice = Lattice_Set(H, W, Sim, Iter, Inlet, PropagationList,
                              CollisionRules, Bin2Dec, FluidCellTypes, 
                              SolidCellTypes, True, SaveParameters, LinkSpeed, 
                              ScalingFactors, GrainSize, ObstacleIndex, Neighborhood_Count,
                              Gas_Sub_Count, Solid_Sub_Count, AdsCoefficient,
                              AdsCapacity, AdsCapacityAccupy, DesCoefficient,
                              DiffCoefficient,RandBankLength,[],TotalSolidWeight)#,SolidWeight,FreeVolume)
        
    elif (choice.lower() == 'l'):
        #----Should be removed at the end------
        print("Loading Previous Data...")
        cf = currentframe()
        filename = getframeinfo(cf).filename
        print("\n\nLoad Simulation should be completed at\n\t",filename,"\n\t Line Number : ",cf.f_lineno)
        #----Should be removed at the end------
        ADD = input('Enter folder name to load data : ') #'out'#############################
        ADD = CurrentDirectory+'/'+ADD
        Name = input('Enter file name to load data : ') # 'Data2000.pkl'##################
        FullPath = ADD+'/'+Name
        Lattice = pickle.load(open(FullPath,"rb"))
        fp = open(ADD+"/Parameters.txt","a+") 
        fp.write("\n\nContinue from file : "+Name)
#        Lattice.Iter.Max = 5000
#        Lattice.Iter.Display = 500
#        Lattice.Iter.Save = 200
        if len(Lattice.Adsorbed[0])<=Lattice.Iter.Max:
            Lattice.Adsorbed = np.append(Lattice.Adsorbed, np.zeros([Lattice.Gas_Sub_Count,Lattice.Iter.Max-len(Lattice.Adsorbed[0])+1]), axis=1)
    elif ((choice == 'C') or (choice == 'c')):
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