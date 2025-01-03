import re
import numpy as np
from viranet.info import metDict, ntpsDict, aaDict, miscDict, ntpsMets, aaMets, N_A, k_atp, k_ppi
from cobra import Model, Reaction, Metabolite
import seaborn as sns
sns.set(context='notebook', style='ticks', font_scale=2)
import matplotlib.pyplot as plt
import pandas as pd

# Parse NCBI gb file to get information
def StringTrick(virusFile, protein):
    # Retrieve translation information from stripped virus file
    tempS1          = [jj for jj, s in enumerate(virusFile) if '                     /gene='+'"'+protein+'"' in s]
    tempSeq      = virusFile[tempS1[0]:]
    tempS2          = [jj for jj, s in enumerate(tempSeq) if '/translation' in s]
    tempSeq      = tempSeq[tempS2[0]:]
    tempS3          = [jj for jj, s in enumerate(tempSeq) if 'gene' in s]
    try:
        tempSeq      = tempSeq[:tempS3[0]]
    except:
        threeUTR = [jj for jj, s in enumerate(tempSeq) if "3'UTR" in s] # HIV has 3'UTR
        tempSeq      = tempSeq[:threeUTR[0]]
    # Clean-up and Store
    regex           = re.compile('[^a-zA-Z]')
    npReg           = re.compile('/translation=')
    protein_sequence  = str(tempSeq)
    protein_sequence  = npReg.sub('',protein_sequence)
    protein_sequence  = regex.sub('',protein_sequence)
    # Cut the redundant information in lower case off
    redundant_start = len(protein_sequence)-1
    for i in protein_sequence:
        if i.islower() == True:
            redundant_start = protein_sequence.index(i)
            break
    protein_sequence = protein_sequence[:redundant_start]
    return protein_sequence
    
# Calculate total AA requirement to generate a viron 
def TotalAACalc(protein_copy_number_dict, virus_file):
    protein_sequence_dict = {}
    # Create protein_sequence_dictionary
    for protein in list(protein_copy_number_dict.keys()):
        protein_sequence_dict[protein] = StringTrick(virus_file, protein)
    # Multiply the protein sequence by its assigned copy number
    totalaa = ''
    for i in range(len(protein_copy_number_dict)):
        protein = list(protein_copy_number_dict.keys())[i]
        totalaa += protein_sequence_dict[protein]*protein_copy_number_dict[protein]
    totalaa_array = np.zeros((20,1))
    for aa in range(len(aaMets)):
        totalaa_array[aa,0] = totalaa.count(aaMets[aa])
    return protein_sequence_dict, totalaa_array

# Actually create cobra reaction object
def CreateCobraReaction(virusFile, Cg, genome_strand, totAA, 
protein_copy_number_dict,protein_sequence_dict, reaction_name):
    for i in metDict.keys():
        if '-L' in metDict[i]:
            metDict[i] = metDict[i].replace('-L', '_L')
        if "[c]" in metDict[i]:
            metDict[i] = metDict[i].replace("[c]", '_c')
	# Define metabolites	
    atp_c = Metabolite(metDict['atp'])
    ctp_c = Metabolite(metDict['ctp'])
    gtp_c = Metabolite(metDict['gtp'])
    utp_c = Metabolite(metDict['utp'])
    ala_c = Metabolite(metDict['A'])
    arg_c = Metabolite(metDict['R'])
    asn_c = Metabolite(metDict['N'])
    asp_c = Metabolite(metDict['D'])
    cys_c = Metabolite(metDict['C'])
    gln_c = Metabolite(metDict['Q'])
    glu_c = Metabolite(metDict['E'])
    gly_c = Metabolite(metDict['G'])
    his_c = Metabolite(metDict['H'])
    ile_c = Metabolite(metDict['I'])
    leu_c = Metabolite(metDict['L'])
    lys_c = Metabolite(metDict['K'])
    met_c = Metabolite(metDict['M'])
    phe_c = Metabolite(metDict['F'])
    pro_c = Metabolite(metDict['P'])
    ser_c = Metabolite(metDict['S'])
    thr_c = Metabolite(metDict['T'])
    trp_c = Metabolite(metDict['W'])
    tyr_c = Metabolite(metDict['Y'])
    val_c = Metabolite(metDict['V'])
    adp_c = Metabolite(metDict['adp'])
    h2o_c = Metabolite(metDict['h2o'])
    h_c   = Metabolite(metDict['h'])
    pi_c  = Metabolite(metDict['Pi'])
    ppi_c = Metabolite(metDict['PPi'])
    # Genome Sequence: initial step is to identify start/end positions in the virus file
    startG  = [jj for jj, s in enumerate(virusFile) if 'ORIGIN' in s]
    endG    = [jj for jj, s in enumerate(virusFile) if '//' in s]
    startG  = int(''.join(map(str,startG)))
    endG    = int(''.join(map(str,endG)))
    startG  = startG + 1
    # Store genome sequence
    regex           = re.compile('[^a-zA-Z]')
    virusGenome     = str(''.join(virusFile[startG:endG]))
    virusGenome     = regex.sub('',virusGenome)
    # [3] Precursor frequency
    # Genome
    # Multiply each NA by its genome copy number
    countA = Cg* virusGenome.count('a')
    countC = Cg* virusGenome.count('c')
    countG = Cg* virusGenome.count('g')
    countU = Cg* virusGenome.count('t')    # Base 'T' is psuedo for base 'U'
    
    if genome_strand==1: # virus whose genome is single-stranded
        antiA = 0
        antiC = 0
        antiG = 0
        antiU = 0
    if genome_strand==2: # virus whose genome is double-stranded
        antiA = countU
        antiC = countG
        antiG = countC
        antiU = countA
    # [4] VBOF Calculations
    # Nucleotides
    # mol.ntps/mol.virus
    V_a = (Cg*(countA + antiA))
    V_c = (Cg*(countC + antiC))
    V_g = (Cg*(countG + antiG))
    V_u = (Cg*(countU + antiU))
    # g.ntps/mol.virus
    G_a = V_a * ntpsDict["atp"]
    G_c = V_c * ntpsDict["ctp"]
    G_g = V_g * ntpsDict["gtp"]
    G_u = V_u * ntpsDict["utp"]
    # Amino Acids
    # mol.aa/mol.virus
    # g.a/mol.virus
    G_aa    = np.zeros((20,1))
    for ii in range(len(aaMets)):
        G_aa[ii,0] = totAA[ii] * aaDict[aaMets[ii]]
    # Total genomic and proteomic molar mass
    M_v = (G_a + G_c + G_g + G_u) + G_aa.sum()
    # Stoichiometric coefficients
    # Nucleotides [mmol.ntps/g.virus]
    S_atp = 1000 * (V_a/M_v)
    S_ctp = 1000 * (V_c/M_v)
    S_gtp = 1000 * (V_g/M_v)
    S_utp = 1000 * (V_u/M_v)
    # Amino acids [mmol.aa/g.virus]
    S_aa    = np.zeros((20,1))
    for ii in range(len(aaMets)):
        S_aa[ii] = 1000 * (totAA[ii]/M_v)
    # Energy requirements
    # Genome: Phosphodiester bond formation products [Pyrophosphate]
    genTemp = (((countA + countC + countG + countU) * k_ppi) - k_ppi)
    genRep  = (((antiA + antiC + antiG + antiU) * k_ppi) - k_ppi)
    genTot  = genTemp + genRep
    V_ppi   = genTot
    S_ppi   = 1000 * (V_ppi/M_v)
    # Protome: Peptide bond formation [ATP + H2O]
    # Note: ATP used in this process is denoated as ATPe/Ae [e = energy version]
    V_Ae = 0
    for i in range(len(protein_copy_number_dict)):
        protein = list(protein_copy_number_dict.keys())[i]
        atp_needed = (len(protein_sequence_dict[protein]))*k_atp - k_atp
        V_Ae+=atp_needed*protein_copy_number_dict[protein]
    S_Ae = 1000 * (V_Ae/M_v)

    # [5] VBOF Reaction formatting and output
    # Left-hand terms: Nucleotides
    # Note: ATP term is a summation of genome and energy requirements
    S_ATP = (S_atp + S_Ae) * -1
    S_CTP = S_ctp * -1
    S_GTP = S_gtp * -1
    S_UTP = S_utp * -1
    # Left-hand terms: Amino Acids
    S_AA    = S_aa * -1
    S_AAf   = dict()
    for ii in range(len(aaMets)):
        S_AAf[aaMets[ii]] = S_AA[ii,0]
    # Left-hand terms: Energy Requirements
    S_H2O = S_Ae * -1
    # Right-hand terms: Energy Requirements
    S_ADP = S_Ae
    S_Pi = S_Ae
    S_H = S_Ae
    S_PPi = S_ppi
    # Create reaction output
    virus_reaction  = Reaction(reaction_name)
    virus_reaction.name = 'Human immunodeficiency virus 1'
    virus_reaction.subsystem = 'Virus Production'
    virus_reaction.lower_bound = 0
    virus_reaction.upper_bound = 1000
    virus_reaction.add_metabolites(({
        atp_c: S_ATP,
        ctp_c: S_CTP,
        gtp_c: S_GTP,
        utp_c: S_UTP,
        ala_c: S_AAf['A'],
        arg_c: S_AAf['R'],
        asn_c: S_AAf['N'],
        asp_c: S_AAf['D'],
        cys_c: S_AAf['C'],
        gln_c: S_AAf['Q'],
        glu_c: S_AAf['E'],
        gly_c: S_AAf['G'],
        his_c: S_AAf['H'],
        ile_c: S_AAf['I'],
        leu_c: S_AAf['L'],
        lys_c: S_AAf['K'],
        met_c: S_AAf['M'],
        phe_c: S_AAf['F'],
        pro_c: S_AAf['P'],
        ser_c: S_AAf['S'],
        thr_c: S_AAf['T'],
        trp_c: S_AAf['W'],
        tyr_c: S_AAf['Y'],
        val_c: S_AAf['V'],
        h2o_c: S_H2O,
        adp_c: S_ADP,
        pi_c:  S_Pi,
        h_c:   S_H,
        ppi_c: S_PPi}))
    
    # Visualize Virus Reaction
    # Amino acids
    plt.figure(dpi=300, figsize=(16,9))
    plt.title('Amino Acid Composition')
    sns.barplot(data=pd.DataFrame({'Estimated Number Per Viron':totAA.flatten(),
              'Amino Acid':aaMets}),
           x='Amino Acid',y='Estimated Number Per Viron', color='#3DA5D9')
    plt.savefig('Virus_production_AA.pdf', dpi=300, bbox_inches='tight')

    # Nucleotides
    plt.figure(dpi=300)
    plt.title('Nucleotide Composition')
    sns.barplot(data={'Nucleotide':['A','C','G','U'],
                  'Estimated Number Per Viron':[countA,countC,countG,countU]},
                x='Nucleotide',y='Estimated Number Per Viron',
           color='#EA7317')
    plt.savefig('Virus_production_NA.pdf', dpi=300, bbox_inches='tight')
    return virus_reaction

