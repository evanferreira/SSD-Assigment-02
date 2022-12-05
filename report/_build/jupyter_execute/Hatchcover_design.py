#!/usr/bin/env python
# coding: utf-8

# # Input data
# Input data is defined below used to run the calculations:
# 
# | Symbol         | Value  | Description                                      |
# | --------------:| ------:|-------------------------------------------------:|
# | L              | 18500  |Length of the hatch cover [mm]                    |
# | B              | 18200  |Breath of hatch cover [mm]                        |
# | W<sub>max</sub>| 200    |Maximum weight of the hatch cover [kg/m²]         |
# | t<sub>c</sub>  | 2.0    |Corrosion adition for plating and stiffeners [mm] |
# | p<sub>w</sub>  | 34.3   |Wave pressure on deck [kN/mm²]                    |
# | S              | 1.1    |Buckling safety factor for hatch covers           |
# | E              | 206000 |Steel Young's Modulus [N/mm²]                     |

# In[1]:


import numpy as np
#Input data from project description and inital definition of the rules
L     = 18500       #Length of hatch cover [mm]
B     = 18200       #Breath of hatch cover [mm]
W_max = 200*L*B/1e6 #Maximum weight of the hatch cover [kg]
tc    = 2.0         #Corrosion adition for plating and stiffeners [mm]
pw    = 34.3        #Wave pressure on deck [kN/mm²]
S     = 1.1         #Buckling safety factor for hatch covers

#Materials yield and allowable stresses
E         = 206e3         #Steel Young's Modulus [N/mm²]
ReH_A     = 235           #Yield stress for mild steel [N/mm²]
ReH_AH32  = 315           #Yield stress for AH32 steel [N/mm²]
ReH_AH36  = 355           #Yield stress for AH36 steel [N/mm²]
sig_a_A   = 0.80*ReH_A    #Allowable normal stress for mild steel [N/mm²]
tal_a_A   = 0.46*ReH_A    #Allowable shear stress for mild steel [N/mm²]
sig_a_32  = 0.80*ReH_AH32 #Allowable normal stress for AH32 steel [N/mm²]
tal_a_32  = 0.46*ReH_AH32 #Allowable shear stress for AH32 steel [N/mm²]
sig_a_36  = 0.80*ReH_AH36 #Allowable normal stress for AH36 steel [N/mm²]
tal_a_36  = 0.46*ReH_AH36 #Allowable shear stress for AH36 steel [N/mm²]

#Materials factor
k_A    = 1.0  #Material factor for mild steel
k_AH32 = 0.78 #Material factor for AH32 steel
k_AH36 = 0.72 #Material factor for AH36 steel


# # Geometry definition

# <div style="text-align: justify"> The cells below calculates the requirements for top plating, stiffeners and primary supporting members (PSM) according to the IACS rules. 
# After checking the minimum scantlings and properties to be used, these structures are defined and used as input for FEM modeling and analysis. On the next section the values will be resumed and presented.
# 
# The cell below calculates the plating requirements and defines its properties. </div>

# In[2]:


#Plating requirements - Considering mild steel
#Net thickness
Mat_pl = "A" #Defining material properties
if Mat_pl == "A":
    ReH_pl   = ReH_A            
    sig_a_pl = sig_a_A 
    tal_a_pl = tal_a_A
    k_pl     = k_A
elif Mat_pl == "AH32":
    ReH_pl   = ReH_AH32            
    sig_a_pl = sig_a_32
    tal_a_pl = tal_a_32
    k_pl     = k_AH32
elif Mat_pl == "AH36":
    ReH_pl   = ReH_AH36            
    sig_a_pl = sig_a_36
    tal_a_pl = tal_a_36
    k_pl     = k_AH36
S_max  = L/3                                    #Spacing of PSM [mm]
Fp     = 1.5                                    #Factor for membrane and bending response
s      = 650e-3                                 #Spacing between stiffeners [m]
FW     = 1.0                                    #Rule coefficient 
tp_min = 15.8*Fp*s*np.sqrt(FW*pw/(0.95*ReH_pl)) #Net thickness of hatch cover top plating [mm]
if tp_min < 10*s or tp_min < 6:
    tp_min = max(10*s,6)

#Plating Net thickness 
tp = tp_min+2.5 #Defining the plating net thickness for FEMAP


# <div style="text-align: justify"> In this step, the ordinary stiffeners requirements are calculated and their final characteristics are set.</div>

# In[3]:


#Ordinary stiffeners requirements
#Ratio hw/tw
Mat_stf = "AH32" #Defining material properties
if Mat_stf == "A":
    ReH_stf   = ReH_A            
    sig_a_stf = sig_a_A 
    tal_a_stf = tal_a_A
    k_stf     = k_A
elif Mat_stf == "AH32":
    ReH_stf   = ReH_AH32            
    sig_a_stf = sig_a_32
    tal_a_stf = tal_a_32
    k_stf     = k_AH32
elif Mat_stf == "AH36":
    ReH_stf   = ReH_AH36            
    sig_a_stf = sig_a_36
    tal_a_stf = tal_a_36
    k_stf     = k_AH36
        
rat_hw_tw  = 15*np.sqrt(235/ReH_stf) #Max. ratio hw/tw for flat bar stiffeners
tw_stf_min = 4                       #Min. web net thickness for stiffeners [mm]

c   = 1.3     #Rule coefficient when plating is stiffened by PSM
ls  = 4625e-3 #Length of the longer side of the plate panel [m]
psi = 1       #Ratio between smallest and largest compressive stress 
ss  = 4550e-3 #Spacing between PSM [m]
m   = c*((1+(ss/ls)**2)**2)*(2.1/(psi+1.1)) #Rule coefficient

#Net section modulus and shear area minimum requirements
w   = (FW*pw*s*ls**2)*1e3/(m*sig_a_stf) #net section modulus of stiffeners [mm³]
Ash = 5*FW*pw*s*ls/tal_a_stf            #net shear area of stiffeners [cm²]

#-------------------------------------------------------------------------------------------------------

#Ordinary stiffeners definition - Angle bar 
hw_stf = 250           #Web height of ordinary stiffener [mm]
bf_stf = 83            #Face plate/flange width of ordinary stiffener [mm]
tw_stf = 7             #Net web thickness of ordinary stiffener [mm]
tf_stf = 11            #Net face plate/flange thickness of ordinary stiffener [mm]
ratio  = hw_stf/tw_stf #Ratio hw/tw for the stiffener

Iw     = ((bf_stf**3) * (hw_stf**2) * 
          (tf_stf * (bf_stf**2 + 2*bf_stf*hw_stf + 4*hw_stf**2) + 3*tw_stf*bf_stf*hw_stf)* 
          1e-6/(12 * (bf_stf + hw_stf)**2))  #Sectorial moment of inertia [cm^6]

Ip     = ((hw_stf**3)*(tw_stf)/3 + 
          hw_stf**2 * bf_stf*tf_stf)*1e-4 #Polar moment of inertia [cm^4]

It     = ((1/3)*(hw_stf*tw_stf**3 + 
                 bf_stf*tf_stf**3 * 
                 (1-0.63*tf_stf/bf_stf)) * 1e-4) #St Venant's moment of inertia [cm^4]

A      =  (hw_stf*tw_stf + 
           bf_stf*tf_stf + 
           tp*(s*1e3))*1e-2 #Cross-sectional area of the stiffener, including face plate [cm²]
l      = ls                 #Stiffener span [m]

stf_na = ((hw_stf*tw_stf)*(hw_stf/2+tp) + 
          (bf_stf*tf_stf)*(tf_stf/2+hw_stf+tp) +
         (tp*s*1e3)*(tp/2))/(A*1e2) #Stiffener neutral axis [mm]

Ia     = ((s*1e3)*tp**3/12 + tw_stf*hw_stf**3/12 + bf_stf*tf_stf**3/12 + 
          (hw_stf*tw_stf)*(tp + hw_stf/2 - stf_na)**2 + 
          (bf_stf*tf_stf)*(tp + hw_stf + tf_stf/2 - stf_na)**2 +
          ((s*1e3)*tp)*(tp/2 - stf_na)**2
         )*1e-4 #Moment of inertia of the stiffener, including face plate [cm^4]
wa     = Ia/(stf_na*1e-1) #Stiffener section modulus [cm³]

#--------------------------------------------------------------
#Stiffener geometric characteristics without the face plate 
A_w      =  (hw_stf*tw_stf + 
             bf_stf*tf_stf)*1e-2 #Cross-sectional area of the stiffener [cm²]

stf_na_w = ((hw_stf*tw_stf)*(hw_stf/2) + 
            (bf_stf*tf_stf)*(tf_stf/2+hw_stf))/(A_w*1e2) #Stiffener neutral axis [mm]

Ia_w     = (tw_stf*hw_stf**3/12 + bf_stf*tf_stf**3/12 + 
            (hw_stf*tw_stf)*(hw_stf/2 - stf_na_w)**2 + 
            (bf_stf*tf_stf)*(hw_stf + tf_stf/2 - stf_na_w)**2)*1e-4 #Moment of inertia of the stiffener [cm^4]
wa_w     = Ia_w/(stf_na_w*1e-1) #Stiffener section modulus [cm³]


# Below it is possible to check the requirements for the PSM. 

# In[4]:


#Primary Supporting Members (PSM) requirements
Mat_psm = "AH32" #Defining material properties
if Mat_psm == "A":
    ReH_psm   = ReH_A            
    sig_a_psm = sig_a_A 
    tal_a_psm = tal_a_A
    k_psm     = k_A
elif Mat_psm == "AH32":
    ReH_psm   = ReH_AH32            
    sig_a_psm = sig_a_32
    tal_a_psm = tal_a_32
    k_psm     = k_AH32
elif Mat_psm == "AH36":
    ReH_psm   = ReH_AH36            
    sig_a_psm = sig_a_36
    tal_a_psm = tal_a_36
    k_psm     = k_AH36
tw_psm_min = 6 #Minimum web net thickness [mm]
tw_psm     = 6 #Defined web net thickness [mm]

#------------------------------------------------------------------------------------------------------

#PSM definition - T bar 
hw        = 1000      #PSM(Primary Supporting Member) depth [mm]
bf_min    = hw*0.4    #(Minimum)Breadth of PSM face plate [mm]
tw_min    = bf_min/30 #(Minimum)Thickness of PSM face plate [mm]
bf        = 400       #Defined breadth of PSM face plate [mm]
tw        = 15        #Defined net thickness of PSM face plate [mm]


# <div style="text-align: justify"> Now, since all the geometry and materials were defined, the cell below calculates the total gross mass (including corrosion addition) of the hatch cover structure.</div>

# In[5]:


#Total mass of the hatch cover
n_PSM = B/(ss*1e3)+1                                    #Number of PSM
n_s   = B/(s*1e3)-n_PSM-1                               #Number of ordinary stiffeners
n_trs = L/(ls*1e3)+1                                    #Number of transversal supports
m_PSM = (((tc+tw)*bf+(tw_psm+tc)*hw)*L*1e-9)*7.85*n_PSM #Mass of PSMs [t]
m_trs = (((tc+tw_psm)*hw+(tw+tc)*bf)*B*1e-9)*7.85*n_trs #Mass of transversal supports
m_pl  = ((tp+tc)*B*L*1e-9)*7.85                         #Mass of top plating [t]
m_stf = (((tw_stf+tc)*hw_stf+
          (tf_stf+tc)*bf_stf)*L*1e-9)*7.85*n_s          #Mass of ordinary stiffeners [t]

M_t   = m_stf+m_pl+m_PSM+m_trs                          #Total mass of the hatch cover [t]


# # Topology of the hatch cover
# This section summarizes the data used and calculated. The cell below shows the hatch cover topology for the design.

# In[6]:


#Summary of inputs
print(f"Maximum spacing between PSM:                  {S_max:.2f} mm")
print(f"Defined spacing between PSM:                  {B/4:.2f} mm")
print(f"Defined height of PSM:                        {hw:.2f} mm")
print(f"Defined spacing between stiffeners:           {s*1000:.2f} mm")
print(f"Number of PSM:                                {n_PSM:.0f}")
print(f"Number of ordinary stiffeners:                {n_s:.0f}")
print(f"Number of transversal supporting members:     {n_trs:.0f}")


# # Minimum requirements - IACS rules
# <div style="text-align: justify"> This section shows the minimum requirements calculated by the code using the IACS rules.</div> 

# In[7]:


#Summary of rule calculations
print(f"Min. net thickness of hatch cover top plating: {tp_min:.2f} mm")
print(f"Max. ratio hw/tw for flat bar stiffeners:      {rat_hw_tw:.2f}")
print(f"Min. net section modulus of stiffeners:        {w:.2f} mm³")
print(f"Min. net shear area of stiffeners:             {Ash:.2f} cm²")
print(f"Min. breadth of PSM face plate:                {bf_min:.2f} mm")
print(f"Min. net thickness of PSM face plate:          {tw_min:.2f} mm")
print(f"Min. net thickness of PSM web:                 {tw_psm_min:.2f} mm")


# # Structure characteristics 
# <div style="text-align: justify"> After calculating the minimum requirements, setting the geometry and material of the structures, below it is shown the resume of the characteristics for top plating, stiffeners and PSM as well as their material properties definition.</div>

# In[8]:


#Summary of hatch cover calculated data
print("-------------------------------------------Top plating------------------------------------------------------")
print(f"Material of the top plating:                  Steel grade {Mat_pl} with yield stress {ReH_pl} MPa")
print(f"Defined net thickness of top plating:         {tp:.2f} mm")
print("-------------------------------------------Stiffeners-------------------------------------------------------")
print(f"Material of the stiffeners:                   Steel grade {Mat_stf} with yield stress {ReH_stf} MPa")
print(f"Web height of stiffener:                      {hw_stf:.2f} mm")
print(f"Face plate/flange width of stiffener:         {bf_stf:.2f} mm")
print(f"Net web thickness of stiffener:               {tw_stf:.2f} mm")
print(f"Net face plate/flange thickness of stiffener: {tf_stf:.2f} mm")
print(f"Sectorial moment of inertia:                  {Iw:.2f} mm\u2076")
print(f"Polar moment of inertia:                      {Ip:.2f} mm\u2074")
print(f"St Venant's moment of inertia:                {It:.2f} mm\u2074")
print(f"Area of the stiffener (including plate):      {A:.2f} cm²")
print(f"Moment of inertia of the stiffener:           {Ia:.2f} cm\u2074")
print(f"Section Modulus of the stiffener:             {wa:.2f} cm³")
print("---------------------------------------------PSM properties-------------------------------------------------")
print(f"Material of the PSM:                          Steel grade {Mat_psm} with yield stress {ReH_psm} MPa")
print(f"Defined breadth of PSM face plate:            {bf:.2f} mm")
print(f"Defined net thickness of PSM face plate:      {tw:.2f} mm")
print(f"Defined net thickness of PSM web:             {tw_psm:.2f} mm")
print("----------------------------------------Hatch cover mass-----------------------------------------------------")
if M_t < W_max/1000:
    print(f"Total gross mass of the hatch cover:          {M_t:.2f} t")
else:
    print(f"Total gross mass of the hatch cover shall be reduced")


# # Finite Element Modeling
# <div style="text-align: justify"> Following LR's procedures and recommendations the hatch cover was modeled. In order to save computational efforts only half of the structure was modeled and boundary conditions were imposed to represent the real physics of the problem. On the next figure the final model with boundary conditions and loads is presented. </div>
# 
# ```{figure} /images/FEM_model.jpg
# ---
# name: FEM_model-fig
# scale: 50
# ---
# ```
# ```{figure} /images/FEM_model_2.jpg
# ---
# name: FEM_model_2-fig
# scale: 45
# ---
# FEM model used on simulation
# ```

# # FEM results
# <div style="text-align: justify"> This section shows the post-processing images and FEM results used to check the buckling criteria. The following measures were checked. 
# <ol>
#     <li>Total deflection of the hatch cover</li>
#     <li>Compressive stress in top plating parallel to stiffeners</li>
#     <li>Compressive stress in top plating perpendicular to stiffeners</li>
#     <li>Compressive stress in stiffeners</li>
#     <li>Shear stress in PSM web panels</li>
#     </ol>
# </div>

# ```{figure} /images/Total_deformation.jpg
# ---
# name: deflection-fig
# ---
# Deflection of the hatch cover
# ```
# 
# ```{figure} /images/Compressive_stress_in_top_plating_parallel_to_stiffeners.jpeg
# ---
# name: stress_in_top_plating_parallel-fig
# ---
# Compressive stress in top plating parallel to stiffeners
# ```
# 
# ```{figure} /images/Compressive_stress_in_top_plating_perpendicular_to_stiffeners.jpeg
# ---
# name: stress_in_top_plating_perpendicular-fig
# ---
# Compressive stress in top plating perpendicular to stiffeners
# ```
# 
# ```{figure} /images/Compressive_stress_in_stiffeners.png
# ---
# name: stress_in_stiffeners-fig
# ---
# Compressive stress in stiffeners
# ```
# 
# ```{figure} /images/Shear_streess_on_PSM_web_panels.jpg
# ---
# name: shear_stress_in_PSM-fig
# ---
# Shear stress in PSM web panels
# ```

# In[9]:


#Input of the results - Enter the compressive stresses
sig_pl_fem_pa  = 103.83 #Stress (parallel to the stiffeners) in the hatch cover plating from FEA [N/mm² (MPa)]
sig_pl_fem_pe  = 108.97 #Stress (perpendicular to stiffeners) in the hatch cover plating from FEA [N/mm² (MPa)]
sig_stf_fem    = 128.93 #Stress in the face plate of the stiffeners from FEA [N/mm² (MPa)]
tal_psm_fem    = 87.476 #Shear stress in the web panels of the PSM from FEA [N/mm² (MPa)]
d_fem          = 63.607 #Total deflection of the PSM from FEA [mm]

#Summary of FEM results
print(f"Stress (parallel to the stiffeners) in the hatch cover plating from FEA:  {sig_pl_fem_pa:.2f} MPa")
print(f"Stress (perpendicular to stiffeners) in the hatch cover plating from FEA: {sig_pl_fem_pe:.2f} MPa")
print(f"Stress in the face plate of the stiffeners from FEA:                      {sig_stf_fem:.2f} MPa")
print(f"Shear stress in the web panels of the PSM from FEA:                       {tal_psm_fem:.2f} MPa")
print(f"Deflection of the PSM from FEA:                                           {d_fem:.2f} mm")


# # Final checks and remarks
# <div style="text-align: justify"> The stresses were check with the criteria and the hatch cover complies with the rule requirements. The materials were chosen aiming to comply with the rules with the lowest grade possible, however steel AH32 was used to comply with some rules requirements, otherwise, using mild steel, the structure would be heavier than the maximum allowed by the input constraints. </div>

# In[10]:


#Maximum deflection
d_max = 0.0056*L #Maximum allowable deflection [mm] 

#------------------------------------------------------------------------------------
#Critical buclking stress check: Compressive stresses in the hatch cover plating
#Stress induced by the bending of PSM, parallel to the direction of stiffeners: 
sig_E1 = 3.6*E*(tp/(s*1e3))**2 #Elastic buckling stress [N/mm²]

#Critical buckling stress [N/mm²]
if sig_E1 > ReH_pl/2:
    sig_C1 = ReH_pl*(1-(ReH_pl/(4*sig_E1)))
else:
    sig_C1 = sig_E1 
sig_c_pa = 0.88*sig_C1/S #Compressive stress threshold

#Stress induced by the bending of PSM, perpendicular to the direction of stiffeners:
sig_E2 = 0.9*m*E*(tp/(s*1e3))**2 #Elastic buckling stress [N/mm²]
#Critical buckling stress [N/mm²]
if sig_E2 > ReH_pl/2:
    sig_C2 = ReH_pl*(1-(ReH_pl/(4*sig_E2)))
else:
    sig_C2 = sig_E2 
sig_c_pe = 0.88*sig_C2/S #Compressive stress threshold

#------------------------------------------------------------------------------------------------

#Critical buclking stress check: Compressive stresses in the face plate of the stiffeners
#Stress induced by the bending of PSM:

sig_E3   = 1e-3*E*Ia/(A*l**2) 
etap     = sig_stf_fem/sig_E1
if 1-etap < 0:
    kp = 0.1
else: 
    kp = 1-etap
C        = ((kp*E*tp**3)*1e-3 /
            (3*s*(1+(1.33*kp*hw_stf*tp**3) /
                  (1000*s*tw_stf**3)))) #Spring stiffness exerted by the top plating
K        = (C*l**4)*1e6/((np.pi**4)*E*Iw)

#Calculating the number of half waves
if K > 0 and K < 4:
    m_w = 1
elif K < 36: 
    m_w = 2
elif K < 144:
    m_w = 3
else: 
    while K >= (m_w**2)*(m_w+1)**2:
        m_w += 1

sig_E4   = ((np.pi**2)*E*Iw/(1e4*Ip*l**2))*(m_w**2+K/m_w**2)+0.385*E*It/Ip

sig_ES   = min(sig_E3, sig_E4) #Elastic buckling stress [N/mm²]
#Critical buckling stress [N/mm²]: 
if sig_ES <= ReH_stf/2: 
    sig_CS = sig_ES 
else: 
    sig_CS = ReH_stf*(1-ReH_stf/(4*sig_ES))

sig_c_stf = 0.88*sig_CS/S #Compressive stress threshold

#-------------------------------------------------------------------------------------

#Critical shear buclking stress check: Shear stresses in the web panels of the PSM:
a  = L*1e-3  #Greater dimension of web panel of PSM [m] 
d  = hw*1e-3 #Smaller dimension of web panel of PSM [m] 
kt = 5.35+4*(a/d)**2

tal_E = 0.9*kt*E*(tw/(1000*d))**2 
if tal_E <= ReH_psm/(2*np.sqrt(3)):
    tal_C = tal_E
else: 
    tal_C = (ReH_psm/np.sqrt(3))*(1-ReH_psm/(4*np.sqrt(3)*tal_E))
    
tal_psm = 0.88*tal_C/S #Sheer stress threshold

#------------------------------------------------------------------------------------------------

#Summary of the results and checks
if d_fem > d_max:
    print(f"Deflection of the hatch cover: {d_fem:.2f} mm \nThreshold: {d_max:.2f} mm - NOK")
else:
    print(f"Deflection of the hatch cover: {d_fem:.2f} mm \nThreshold: {d_max:.2f} mm - OK")
    
if sig_pl_fem_pa > sig_c_pa:
    print(f"Stress (parallel to stiffeners) in the hatch cover plating: {sig_pl_fem_pa:.2f} MPa \nThreshold: {sig_c_pa:.2f} MPa - NOK")
else:
    print(f"Stress (parallel to stiffeners) in the hatch cover plating: {sig_pl_fem_pa:.2f} MPa \nThreshold: {sig_c_pa:.2f} MPa - OK")
        
if sig_pl_fem_pe > sig_c_pe:
    print(f"Stress (perpendicular to stiffeners) in the hatch cover plating {sig_pl_fem_pe:.2f} MPa \nThreshold: {sig_c_pe:.2f} MPa - NOK")
else:
    print(f"Stress (perpendicular to stiffeners) in the hatch cover plating {sig_pl_fem_pe:.2f} MPa \nThreshold: {sig_c_pe:.2f} MPa - OK")
    
if sig_stf_fem > sig_c_stf:
    print(f"Stress in stiffeners face plate: {sig_stf_fem:.2f} MPa \nThreshold: {sig_c_stf:.2f} MPa - NOK")
else:
    print(f"Stress in stiffeners face plate: {sig_stf_fem:.2f} MPa \nThreshold: {sig_c_stf:.2f} MPa - OK")

if tal_psm_fem > tal_psm:
    print(f"Shear stress in the PSM web panels: {tal_psm_fem:.2f} MPa \nThreshold: {tal_psm:.2f} MPa - NOK")
else:
    print(f"Shear stress in the PSM web panels: {tal_psm_fem:.2f} MPa \nThreshold: {tal_psm:.2f} MPa - OK")

