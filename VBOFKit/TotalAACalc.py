from VBOFKit import StringTrick
import numpy as np
from viranet.info import metDict, ntpsDict, aaDict, miscDict, ntpsMets, aaMets, N_A, k_atp, k_ppi
# Total AA calculation from the aggregate protein sequence object created by Stringtrick function
def TotalAACalc(protein_copy_number_dict, virus_file):
    protein_sequence_dict = {}
    # Create protein_sequence_dictionary
    for protein in list(protein_copy_number_dict.keys()):
        protein_sequence_dict[protein] = StringTrick.StringTrick(virus_file, protein)
    # Multiply the protein sequence by its assigned copy number
    totalaa = ''
    for i in range(len(protein_copy_number_dict)):
        protein = list(protein_copy_number_dict.keys())[i]
        totalaa += protein_sequence_dict[protein]*protein_copy_number_dict[protein]
    totalaa_array = np.zeros((20,1))
    for aa in range(len(aaMets)):
        totalaa_array[aa,0] = totalaa.count(aaMets[aa])
    return protein_sequence_dict, totalaa_array
    

