from Initialize import Initialize
import numpy as np

print ("\n"*80)

StartChoice = 'N'
MH = 4e-5 # Macroscopic Value of the H
H = 21
MW = MH   # Macroscopic Value of the W
Gas_Sub_Count=1########2############
Flow = np.array([3])############[[20],[10]]############
Solid_Sub_Count=1########2############
SolidCount = np.array([10])##########np.array([20,5])######
IsoThermal=True
GrainSize = 1 ##################

MT = 1.2e-5  # Macroscopic Value of the Final Simulation Time
SaveIter = 1000 # Save iteration number
DispIter = 100 # Display iteration number

MGPC = sum(Flow)*7 # Maximum Gas Particles in each cell
GMM = 0.018 # Gas Molar Mass (Kg/mol)
GD = 0.58 # Gas Adsorbate Density (Kg/m^3)
MZ = 1.6e-7 # Macroscopic Thickness of the Cell (m)

M_AdsCoefficient = np.array([[1.47e7]]) # Adsorption rate constant (m3 mol-1 s-1) [Solid_Sub_Count,Gas_Sub_Count]
M_DesCoefficient = np.array([[2.94e6]]) # Desorption rate constant (s-1) [Solid_Sub_Count,Gas_Sub_Count]
M_DiffCoefficient = np.array([[3.92e-8]]) # Solid diffusion coefficient (m^2/s) Solid_Sub_Count,Gas_Sub_Count]

M_AdsCapacity = np.array([0.32]) #Saturation adsorption Weight (Kg/Kg) [Solid_Sub_Count,Gas_Sub_Count]
SD = np.array([670]) # (Kg/m3) Solid Density
GD = np.array([0.58]) # (Kg/m3) Gas Density
AdsCapacityAccupy = np.array([1]) #[Gas_Sub_Count]

Energy = []

ADD = 'out' # Folder name to save data
Name = 'Data' # File name to save data

SimParams = [MH, H, MW, Gas_Sub_Count, Flow, Solid_Sub_Count, SolidCount, 
             IsoThermal, GrainSize, MT, SaveIter, DispIter, GMM, GD, MZ,
             M_AdsCoefficient, M_DesCoefficient, M_DiffCoefficient, 
             M_AdsCapacity, SD, GD, Energy, ADD, Name]

ObstacleChoice = 2
#2 # SquarePorousObstacle
#"0: No Obstacle"
#"1: Load from previous simulation file"
#"2: SquarePorousObstacle"
#"3: SingleRoundObstacle"
#"4: SingleSquareObstacle"
#"5: QSGS"

Case = 'S1'

i1 = 1
j1 = 1
Height = H
Width = H
if Case == 'S1':
    NumColumns = 4
    BlockSize = int(H/20)
elif Case == 'S2':
    NumColumns = 8
    BlockSize = int(H/20)
elif Case == 'S3':
    NumColumns = 12
    BlockSize = int(H/20)
elif Case == 'S4':
    NumColumns = 8
    BlockSize = int(3*H/40)
elif Case == 'S5':
    NumColumns = 4
    BlockSize = int(3*H/20)
else:
    NumColumns = 4 #################
    BlockSize = 1#int(3*H/20) #################

    
BlockDistance = int((H-BlockSize*NumColumns)/(NumColumns+1))
Shift = [int(1.5*BlockDistance),int(0.5*BlockDistance)]
            
                
ObstacleParams = [i1, j1, Height, Width, NumColumns, BlockSize, BlockDistance, Shift, Case]
            
Lattice = Initialize(StartChoice, SimParams, ObstacleChoice, ObstacleParams)
Lattice.Simulate()
Lattice.Output(3)