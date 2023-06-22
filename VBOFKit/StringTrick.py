import re
import numpy as np
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
