# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 19:27:04 2025

@author: zacha
"""

from biocrnpyler import *
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

def run_sim(timepoints, x0, CRN):
    R = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, 
                                             stochastic=True)
    return R['protein_CFP_degtagged'].iloc[-1]

def rna_se_biocrnpyler(ktl, offset, color):
    
    # Parameters (found in ctRSD_simulator_210.py)
    params = {"ktx":0.013, "kdil":0.00075, "kcat":0.25/60} # kdil not from source, use what value? 
    #should there be a different kdil for dsRNA
    
    krsd = 1e3/1e9
    krev = 270/1e9
    ktl = ktl # not from source, use what value?
    
    # Species
    input_rna = Species('input_rna', material_type='rna')
    output_rna = Species('output_rna', material_type='rna')
    
    ugate = Species('ugate', material_type='rna')
    gate = Species('gate', material_type='rna')
    input_gate = Species('input_gate', material_type='rna')
    
    CFP = Species('CFP', material_type='protein', attributes=['degtagged'])
    
    ribozyme = Enzyme('ribozyme', substrates = [ugate], products = [gate]) #not a protein but defaults to protein material_type
    
    # DNAassemblies and Mixture 
    assembly_ugate = DNAassembly(name='ugate', promoter='strong', rbs=None,
                                       transcript=ugate)
    assembly_input = DNAassembly(name='input', promoter='strong', rbs=None, 
                                 transcript=input_rna)
    
    dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)
    
    global_mechanisms = {"dilution":dilution_mechanism}
    
    M = SimpleTxTlExtract(name="txtl", parameters = params, species=[CFP, input_gate, output_rna], 
                          global_mechanisms=global_mechanisms,
                          components=[assembly_ugate, assembly_input, ribozyme])
    
    # Strand Displacement and Translation
    SD = Reaction.from_massaction([input_rna, gate], [input_gate, output_rna], k_forward=krsd, 
                                  k_reverse=krev)
    Tl = Reaction.from_massaction([input_gate], [CFP, input_gate], k_forward=ktl)
    
    # File Creation and Simulation
    CRN = M.compile_crn()
    CRN.add_reactions([SD, Tl])
    CRN.write_sbml_file(f'rna-strand-exchange.xml') #saving CRN as sbml
    
    with open('temp_CRN_EQNs.txt', 'w') as f:
        f.write(CRN.pretty_print(show_rates = True, show_keys = True))
    
    print('CRN Compiled')
    
    timepoints = np.linspace(0, 12000, 500)
    
    x0 = {assembly_ugate.dna:25, assembly_input.dna: 25, 'protein_ribozyme':1}
    
    data = []
    
    for i in range(10000):
        R =  CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, 
                                                 stochastic=True)
        data.append(R['protein_CFP_degtagged'].iloc[-1])
    
    y, x = np.histogram(data, bins=200)
    x_mid = 0.5 * (x[1:] + x[:-1])
    
    plt.fill_between(x_mid, offset, y + offset, color=color, alpha=0.3)
    plt.plot(x_mid, y + offset, color=color, lw=1.2)
    plt.plot([1, 2e4], [offset, offset], color='k')

#parameter list for function
ktl_lst = [0.075, 0.05, 0.015, 0.0075, 0.0025, 0.001]
color_lst = ['#123524', '#1D5638', '#277A4C', '#32A160', '#4DCB87', '#7FEFB2']
offset_lst = np.arange(len(ktl_lst)) * 200

#Initialize Figure
plt.figure()

for j in range(len(ktl_lst)):
    rna_se_biocrnpyler(ktl_lst[j], offset_lst[j], color_lst[j])
    
plt.title('CFP Stochastic Distribution')
plt.xscale('log')
plt.yticks([])
plt.xlabel('CFP Signal')
plt.ylabel('Cell Count')

plt.show()
