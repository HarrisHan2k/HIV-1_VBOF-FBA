import cobra
import pandas as pd
# Optimization
def Optimize(HVM,HostRxn, VirusRxn, solver):
    # Optimize the model
    hostIdx = HostRxn
    virusIdx    = VirusRxn  # usually the virus-generation reaction is added at the end the model  
    objIdx      = [hostIdx,virusIdx]
    # Set the solver of HVM
    HVM.solver = solver
    # [2] State Optimisations
    # Identify the reactions
    hostObj     = HVM.reactions[hostIdx]
    virusObj    = HVM.reactions[virusIdx]
    # Record the bounds
    hostLb      = HVM.reactions[hostIdx].bounds[0]
    hostUb      = HVM.reactions[hostIdx].bounds[1]
    virusLb     = HVM.reactions[virusIdx].bounds[0]
    virusUb     = HVM.reactions[virusIdx].bounds[1]
    # Zero-bound virus reaction
    HVM.reactions[virusIdx].bounds = (0,0)
    # Optimize host fuction
    HVM.objective = hostObj
    hostSol     = HVM.optimize(objective_sense='maximize')
    hostX       = hostSol.fluxes                 # Host-optimal flux vector
    hostF       = hostSol.objective_value                 # Optima value for host objective
    # Return virus reaction to correct bounds
    HVM.reactions[virusIdx].bounds = (virusLb, virusUb)
    # Virus optimization
    # Zero-bound host reaction
    HVM.reactions[hostIdx].bounds = (0,0)
    # Record virus optima
    HVM.objective = virusObj
    virusSol    = HVM.optimize(objective_sense='maximize')
    virusX      = virusSol.fluxes                # Virus-optimal flux vector
    virusF      = virusSol.objective_value                # Optima value for virus objective
    # Return host reaction to correct bounds
    HVM.reactions[hostIdx].bounds = (hostLb,hostUb)
    return (hostF, hostX, virusF, virusX)


def Knockout(HVM,HostRxn, VirusRxn, solver):
    import numpy as np
    import cobra
    "Reaction Knockouts"
    # [1] Initial Setup
    # Identify the host objective reaction
    try:
        intTest = int(HostRxn)
        hostIdx = HostRxn
    except:
        print('HostRxn intTest failed')
        for ii in range(len(HVM.reactions)):
            if HostRxn in str(HVM.reactions[ii]):
                hostIdx = ii
    # Ensure virus is objective
    virusIdx    = VirusRxn
    # Set the solver of HVM
    HVM.solver = solver
    virusObj    = HVM.reactions[virusIdx]
    HVM.objective = virusObj
    # Variable creation to hold virus optima
    koVirus     = np.zeros((len(HVM.reactions),1))
    # [2] Knockout Analysis
    # Record the host bounds
    hostLb      = HVM.reactions[hostIdx].lower_bound
    hostUb      = HVM.reactions[hostIdx].upper_bound
    # Initiate loop
    HVM_ko = HVM
    HVM_ko.reactions[hostIdx].lower_bound = 0
    HVM_ko.reactions[hostIdx].upper_bound = 0
    HVM_ko.optimize(objective_sense='maximize')
    # Ensure the host objective is set to zero bounds
    # Single_gene_ko for virus optima 
    single_gene_ko = cobra.flux_analysis.single_gene_deletion(HVM_ko)
    # Single_reaction_ko for virus optima
    single_reaction_ko = cobra.flux_analysis.single_reaction_deletion(HVM_ko)
    # Change reaction items to string type and remove brackets
    single_reaction_ko['ids'] = single_reaction_ko['ids'].apply(str)
    def RemoveBrackets(i):
        i = i[2:len(i)-2]
        return i
    single_reaction_ko['ids'] = single_reaction_ko['ids'].apply(RemoveBrackets)
    # [3] Outputs
    return single_gene_ko, single_reaction_ko

def AnnotateGene(cheatsheet, KOresults):
    gene_symbol_list = []
    for gene_id in KOresults['ids']:
        gene_symbol = cheatsheet[cheatsheet[cheatsheet.columns[0]]==int(list(gene_id)[0])].values
        if len(gene_symbol)==0:
            gene_symbol_list.append('Unmatched')
            continue
        else:
            gene_symbol_list.append(gene_symbol[0,1])
    KOresults['gene symbol'] = gene_symbol_list
    
def Compare(HVM_remove_obj, objIdx, hostX,virusX):
    import numpy as np
    import pandas as pd
    "Host-Virus Comparison"
    # [1] Initial Setup
    # Numpy Conversion
    hostX   = np.array(hostX)
    virusX  = np.array(virusX)
    # Remove the objective reactions from both flux vectors
    hostXd  = np.delete(hostX,objIdx)
    virusXd = np.delete(virusX,objIdx)
    # Convert flux vectors to absolute and normalise to summation of vector
    pHOS    = (hostXd / sum(np.absolute(hostX))) * 100
    pVOS    = (virusXd / sum(np.absolute(virusX))) * 100
    # Convert to suitable numpy array
    pHOS    = np.array(pHOS)
    pVOS    = np.array(pVOS)
    # [2] Reaction Statistics
    # Variables for calculations
    tol     = 1.05                                          # Regulated tolerance
    e       = 1e-06                                         # Threshold for 'on'
    ne      = e * -1                                        # Negative threshold
    inf     = np.inf
    # Reaction states: Upregulated; Downregulated; Activated; Inactivated; Reversed
    regulation = []
    # Initiate loop
    for ii in range(len(pHOS)):
        # Upregulated
        if (pVOS[ii] > e) and (pHOS[ii] > e) and ((pVOS[ii] / pHOS[ii]) > tol) and ((pVOS[ii] / pHOS[ii]) < inf):
            regulation.append('Upregulated')
            continue
        elif (pVOS[ii] < ne) and (pHOS[ii] < ne) and ((pVOS[ii] / pHOS[ii]) > tol) and ((pVOS[ii] / pHOS[ii]) < inf):
            regulation.append('Downregulated')
            continue
        # Downregulated
        if (pVOS[ii] > e) and (pHOS[ii] > e) and ((pHOS[ii] / pVOS[ii]) > tol) and ((pHOS[ii] / pVOS[ii]) < inf):
            regulation.append('Downregulated')
            continue
        elif (pVOS[ii] < ne) and (pHOS[ii] < ne) and ((pHOS[ii] / pVOS[ii]) > tol) and ((pHOS[ii] / pVOS[ii]) < inf):
            regulation.append('Upregulated')
            continue
        # Activated
        if (np.absolute(pVOS[ii]) > e) and (np.absolute(pHOS[ii]) < e):
            regulation.append('Activated')
            continue
        # Inactivated
        if (np.absolute(pVOS[ii]) < e) and (np.absolute(pHOS[ii]) > e):
            regulation.append('Inactivated')
            continue
        # Reversed
        if (pVOS[ii] > e and pHOS[ii] < ne) or (pVOS[ii] < ne and pHOS[ii] > e):
            regulation.append('Reversed')
            continue
        regulation.append('No difference')
    # Outputs
        # hvmComp
    hvmComp     = np.vstack((pHOS,pVOS))
    hvmComp     = hvmComp.transpose()
        # Regulation
        # Get reaction names
    reactions = []
    for i in HVM_remove_obj.reactions:
        reactions.append(i.id)
    regulation = pd.DataFrame({'Reaction':reactions,
                               'Regulation':regulation})
    hvmComp = pd.DataFrame(hvmComp,
    	index=reactions, columns=['Host optima', 'Virus Optima'])
    return(hvmComp, regulation, hostXd, virusXd)

def FVAEnforce(HVM,virusIdx,hostIdx,fraction_of_optimum, virus_optimum):
    hostObj     = HVM.reactions[hostIdx]
    # Record the bounds
    hostLb      = HVM.reactions[hostIdx].bounds[0]
    hostUb      = HVM.reactions[hostIdx].bounds[1]
    virusLb     = HVM.reactions[virusIdx].bounds[0]
    virusUb     = HVM.reactions[virusIdx].bounds[1]
    # Host optimization
    HVM.objective = hostIdx
    HVM.reactions[virusIdx].bounds = (0,0)
    # Perform flux varability analysis
    varHost     = cobra.flux_analysis.flux_variability_analysis(HVM,loopless=True, fraction_of_optimum=fraction_of_optimum)
    # Return virus reaction to correct bounds
    HVM.reactions[virusIdx].bounds = (virusLb, virusUb)
    # Enforcement loop
    enfVirus = []
    # Ensure the host objective is set to zero bounds set virus reaction as objective
    HVM.objective = virusIdx
    HVM.reactions[hostIdx].bounds = (0,0)
    for ii in range(len(HVM.reactions)):
        # Store the bounds of the reaction [ii]
        rxnLb = HVM.reactions[ii].bounds[0]
        rxnUb = HVM.reactions[ii].bounds[1]
        # Alter reacton [ii] bounds to match median host-derived flux from FVA: varHost[ii]
        maxFVA      = varHost.loc[HVM.reactions[ii].id]['maximum']
        minFVA      = varHost.loc[HVM.reactions[ii].id]['minimum']
        medianFVA   = maxFVA - ((maxFVA - minFVA) / 2)
        HVM.reactions[ii].bounds = (medianFVA,medianFVA)
        # Record and store the virus optima
        vSol = HVM.optimize(objective_sense='maximize')
        enfVirus.append(vSol.objective_value)
        # Return reaction [ii] to it's original bounds
        HVM.reactions[ii].bounds = (rxnLb, rxnUb)
    # Return the host objective
    HVM.reactions[hostIdx].bounds = (hostLb, hostUb)
    reactions = []
    for i in HVM.reactions:
        reactions.append(i.id)
    enfVirus = pd.DataFrame({'Reaction':reactions,
                              'Enforced growth':enfVirus})
    enfVirus['Enforce efficiency'] = enfVirus['Enforced growth']/virus_optimum
    return enfVirus
