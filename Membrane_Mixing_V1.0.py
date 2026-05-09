# =============================================================================
# TITLE:      Membrane Separation Model (Perfect Mixing)
# MODULE:     DWSIM Python Custom Model
# AUTHOR:     ChemPhys Modeling Lab
# VERSION:    1.0.0
# DATE:       May 2026
# CONTACT:    chemphys.modeling@gmail.com
# WEBSITE:    https://chemphysmodeling.gumroad.com/l/Membrane_Model_PerfectMixing
#
# DESCRIPTION:
# This script calculates the permeate and retentate flow rates and 
# compositions for gas separation based on the Perfect Mixing assumption.
# Developed for use within the DWSIM process simulator environment.
#
# STATUS:
# This is an initial functional release. Code refactoring and
# comprehensive commenting are planned for the next version.
# =============================================================================

def calculate_unit_output(P_in_Pa, P_out_Pa, Area_cm2, compMolFlow_in_mols, Permeance_GPU) :
    # Library
    import numpy as np
    from scipy.optimize import least_squares

    def membrane_mixing(vars, P_in_Pa, P_out_Pa, compMolFlow_in_mols, Area_cm2, Permeance_GPU):
    
        permMolFlow_mols = {}
        PPperm_cmHg = {}
        PPret_cmHg = {}

        i = 0
        for compName in compMolFlow_in_mols.keys():
        
            permMolFlow_mols[compName] = vars[i]
            PPret_cmHg[compName] = vars[i+1]
            PPperm_cmHg[compName] = vars[i+2]

            i = i + 3

        totalMolFlow_mols = sum(compMolFlow_in_mols.values())
        totalPermMolFlow_mols = sum(permMolFlow_mols.values()) + totalMolFlow_mols * 1e-9
        totalRetMolFlow_mols = totalMolFlow_mols - totalPermMolFlow_mols + totalMolFlow_mols * 1e-9

        MolFrac_in = {}
        for compName in compMolFlow_in_mols.keys():
            MolFrac_in[compName] = compMolFlow_in_mols[compName] / totalMolFlow_mols
    
        residuals = []
        for compName in compMolFlow_in_mols.keys():
            if Permeance_GPU[compName] == 0 or compMolFlow_in_mols[compName] == 0:
                eq1 = permMolFlow_mols[compName]
                eq2 = PPret_cmHg[compName] - compMolFlow_in_mols[compName] / (compMolFlow_in_mols[compName]+1e-6)
                eq3 = PPperm_cmHg[compName]
            else : 
                eq1 = (- permMolFlow_mols[compName]*22400 + Permeance_GPU[compName]* 1e-6 * (PPret_cmHg[compName] - PPperm_cmHg[compName]) * Area_cm2) / (compMolFlow_in_mols[compName] /100)
                eq2 = (- PPret_cmHg[compName] + P_in_Pa * 0.00075 * (compMolFlow_in_mols[compName] - permMolFlow_mols[compName])  / totalRetMolFlow_mols ) / (P_in_Pa*0.00075*MolFrac_in[compName])
                eq3 = (- PPperm_cmHg[compName] + P_out_Pa * 0.00075 * permMolFlow_mols[compName] / totalPermMolFlow_mols)/(P_out_Pa*0.00075)
        
            residuals.append(eq1)
            residuals.append(eq2)
            residuals.append(eq3)
        
        return residuals

    N_comp = len(list(compMolFlow_in_mols.keys()))

    TotalMolFlow_mols = 0
    PPin_cmHg = {}
    TotalMolFlow_mols = sum(compMolFlow_in_mols.values())

    for compName in compMolFlow_in_mols.keys():
        PPin_cmHg[compName] = P_in_Pa * 0.00075 * compMolFlow_in_mols[compName] / TotalMolFlow_mols

    permMolFlow_mols2 = {}
    PPret_cmHg2 = {}
    PPperm_cmHg2 = {}

    PPret_cmHg3 = {}
    PPperm_cmHg3 = {}

    Number_comps = len(compMolFlow_in_mols.keys())
    
    initial_guess_mixed = []
    for compName in compMolFlow_in_mols.keys():
        initial_guess_mixed.append(compMolFlow_in_mols[compName]/100)
        initial_guess_mixed.append(P_in_Pa * 0.00075 * compMolFlow_in_mols[compName] / TotalMolFlow_mols)
        initial_guess_mixed.append(P_out_Pa * 0.00075 * Permeance_GPU[compName] / sum(Permeance_GPU.values()))
    
    bounds_mixed = []
    lower_bounds = []
    upper_bounds = []

    for compName in compMolFlow_in_mols.keys():
        lower_bounds.append(0.0)
        upper_bounds.append(compMolFlow_in_mols[compName])
        lower_bounds.append(0.0)
        upper_bounds.append(P_in_Pa * 0.00075+1e-11)
        lower_bounds.append(0.0)
        upper_bounds.append(P_out_Pa * 0.00075+1e-11)

    bounds_mixed = (lower_bounds, upper_bounds)

    results_mixed = least_squares(membrane_mixing, initial_guess_mixed, bounds=bounds_mixed,args=(P_in_Pa, P_out_Pa,compMolFlow_in_mols, Area_cm2, Permeance_GPU))

    i = 0
    for compName in compMolFlow_in_mols.keys():
        permMolFlow_mols2[compName] = results_mixed.x[i]
        PPret_cmHg2[compName] = results_mixed.x[i+1]
        PPperm_cmHg2[compName] = results_mixed.x[i+2]
                
        i = i + 3
          
    compMolFlow_Ret_mols = {}
    compMolFlow_Perm_mols = {}
    for compName in compMolFlow_in_mols.keys():
        compMolFlow_Perm_mols[compName] = permMolFlow_mols2[compName]
        compMolFlow_Ret_mols[compName] = compMolFlow_in_mols[compName] - compMolFlow_Perm_mols[compName]
    
    test = 0

    return compMolFlow_Ret_mols, compMolFlow_Perm_mols

def get_inletData_Mock() :

    P_in_Pa = 200000 
    P_out_Pa = 1000
    Area_cm2 = 10*100*100 #[kJ/kg]
    Temperature = 300
    N_cell = 50
        
    compMolFlow_in_mols = {}
    compMolFlow_in_mols["Carbon dioxide"] = 0.2 #[mol/s]
    compMolFlow_in_mols["Nitrogen"] = 0.79 #[mol/s]
    compMolFlow_in_mols["Hydrogen"] = 0.01 #[mol/s]

    Permeance_GPU = {}
    Permeance_GPU["Carbon dioxide"] = 1000 #[mol/s]
    Permeance_GPU["Nitrogen"] = 1 #[mol/s]
    Permeance_GPU["Hydrogen"] = 0 #[mol/s]

    return  P_in_Pa, P_out_Pa, Area_cm2, Temperature, compMolFlow_in_mols, Permeance_GPU, N_cell

if __name__ == "__main__":

    P_in_Pa, P_out_Pa, Area_cm2, Temperature, compMolFlow_in_mols,Permeance_GPU, N_cell = get_inletData_Mock()

    compMolFlow_Ret_mols, compMolFlow_Perm_mols = calculate_unit_output(P_in_Pa, P_out_Pa, Area_cm2, compMolFlow_in_mols, Permeance_GPU)

else :
    import System
    from DWSIM.Thermodynamics import *

    inlet = ims1
    Flowsheet.WriteMessage('import start')
    P_in_Pa = inlet.Phases[0].Properties.pressure
    Temperature = inlet.Phases[0].Properties.temperature

    Area_cm2 = Area_m2 * 10000
    P_out_Pa = 1000 * Pressure_out_kPa

    compMolFlow_in_mols = {}
    for comp in inlet.Phases[0].Compounds.Values:
        compMolFlow_in_mols[comp.Name]  = comp.MolarFlow
    
    Permeance_GPU = {}
    for compName in compMolFlow_in_mols.keys():
        variablename = str('GPU_' + compName)
        variablename = variablename.replace(' ', '_')
        temp = eval(variablename)
        if temp == None:
            Flowsheet.WriteMessage(variablename + 'is null')
            temp = 1e-9
        if temp == 0:
            Flowsheet.WriteMessage(variablename + 'is 0, set to 1e-9')
            temp = 1e-9
        Permeance_GPU[compName] = temp
        Flowsheet.WriteMessage(str(Permeance_GPU[compName]))
        
    compMolFlow_Ret_mols, compMolFlow_Perm_mols = calculate_unit_output(P_in_Pa, P_out_Pa, Area_cm2, compMolFlow_in_mols, Permeance_GPU)

    Flowsheet.WriteMessage('end calcoutput')
    
    outlet_ret = oms1
    outlet_perm = oms2
    
    outlet_ret.Clear()
    outlet_ret.ClearAllProps()
    outlet_perm.Clear()
    outlet_perm.ClearAllProps()

    outlet_ret.Phases[0].Properties.pressure = P_in_Pa
    outlet_ret.Phases[0].Properties.temperature = Temperature
    outlet_perm.Phases[0].Properties.pressure = P_out_Pa
    outlet_perm.Phases[0].Properties.temperature = Temperature

    ret_totalflow = 0
    for comp in inlet.Phases[0].Compounds.Values:
        ret_totalflow += compMolFlow_Ret_mols[comp.Name]
    outlet_ret.Phases[0].Properties.molarflow = ret_totalflow

    perm_totalflow = 0
    for comp in inlet.Phases[0].Compounds.Values:
        perm_totalflow += compMolFlow_Perm_mols[comp.Name]
    outlet_perm.Phases[0].Properties.molarflow = perm_totalflow

    for comp in outlet_ret.Phases[0].Compounds.Values:
        comp.MoleFraction = compMolFlow_Ret_mols[comp.Name]/ret_totalflow
    for comp in outlet_perm.Phases[0].Compounds.Values:
        comp.MoleFraction = compMolFlow_Perm_mols[comp.Name]/perm_totalflow

 
   

