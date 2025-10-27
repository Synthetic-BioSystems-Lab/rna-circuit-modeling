# -*- coding: utf-8 -*-
"""
Created on Wed Oct  22 17:19:04 2025

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
    params = {"ktx":0.013, "kdil":0.00075, "kcat":0.25/60, 
              ParameterKey(mechanism=None, part_id='rna_ugate', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_gate', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_input_gate', name='kdil'):0.000375
              } # kdil not from source, value? 
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
    CRN.write_sbml_file(f'ctRSD.xml') #saving CRN as sbml
    
    with open('temp_CRN_EQNs.txt', 'w') as f:
        f.write(CRN.pretty_print(show_rates = True, show_keys = True))
    
    print('CRN Compiled')
    
    timepoints = np.linspace(0, 12000, 500)
    
    x0 = {assembly_ugate.dna:24, assembly_input.dna: 24, 'protein_ribozyme':1}
    
    data = []
    
    for i in range(100):
        R =  CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, 
                                                 stochastic=True)
        data.append(R['protein_CFP_degtagged'].iloc[-1])
    
    y, x = np.histogram(data, bins=200)
    x_mid = 0.5 * (x[1:] + x[:-1])
    
    plt.fill_between(x_mid, offset, y + offset, color=color, alpha=0.3)
    plt.plot(x_mid, y + offset, color=color, lw=1.2)
    plt.plot([1, 2e4], [offset, offset], color='k')

def rna_se_biocrnpyler_half_and_half(ktl_0, ktl_1, offset, color):
    
    # Parameters (found in ctRSD_simulator_210.py)
    params = {"ktx":0.013, "kdil":0.00075, "kcat":0.25/60, 
              ParameterKey(mechanism=None, part_id='rna_ugate_0', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_gate_0', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_input_gate_0', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_ugate_1', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_gate_1', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_input_gate_1', name='kdil'):0.000375
              } # kdil not from source, value? 
    #should there be a different kdil for dsRNA
    
    krsd = 1e3/1e9
    krev = 270/1e9
    ktl_0 = ktl_0 # not from source, use what value?
    ktl_1 = ktl_1
    
    # Species
    input_rna = Species('input_rna', material_type='rna')
    output_rna = Species('output_rna', material_type='rna')
    
    ugate_0 = Species('ugate_0', material_type='rna')
    ugate_1 = Species('ugate_1', material_type='rna')
    gate_0 = Species('gate_0', material_type='rna')
    gate_1 = Species('gate_1', material_type='rna')
    input_gate_0 = Species('input_gate_0', material_type='rna')
    input_gate_1 = Species('input_gate_1', material_type='rna')
    
    CFP = Species('CFP', material_type='protein', attributes=['degtagged'])
    
    ribozyme_0 = Enzyme('ribozyme_0', substrates = [ugate_0], products = [gate_0]) #not a protein but defaults to protein material_type
    ribozyme_1 = Enzyme('ribozyme_1', substrates = [ugate_1], products = [gate_1])
    
    # DNAassemblies and Mixture 
    assembly_ugate_0 = DNAassembly(name='ugate_0', promoter='strong', rbs=None,
                                       transcript=ugate_0)
    assembly_ugate_1 = DNAassembly(name='ugate_1', promoter='strong', rbs=None,
                                       transcript=ugate_1)
    assembly_input = DNAassembly(name='input', promoter='strong', rbs=None, 
                                 transcript=input_rna)
    
    dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)
    
    global_mechanisms = {"dilution":dilution_mechanism}
    
    M = SimpleTxTlExtract(name="txtl", parameters = params, species=[CFP, input_gate_0, input_gate_1, output_rna], 
                          global_mechanisms=global_mechanisms,
                          components=[assembly_ugate_0, assembly_ugate_1, assembly_input, 
                                      ribozyme_0, ribozyme_1])
    
    # Strand Displacement and Translation
    SD_0 = Reaction.from_massaction([input_rna, gate_0], [input_gate_0, output_rna], k_forward=krsd, 
                                  k_reverse=krev)
    SD_1 = Reaction.from_massaction([input_rna, gate_1], [input_gate_1, output_rna], k_forward=krsd, 
                                  k_reverse=krev)
    Tl_0 = Reaction.from_massaction([input_gate_0], [CFP, input_gate_0], k_forward=ktl_0)
    Tl_1 = Reaction.from_massaction([input_gate_1], [CFP, input_gate_1], k_forward=ktl_1)
    
    # File Creation and Simulation
    CRN = M.compile_crn()
    CRN.add_reactions([SD_0, Tl_0, SD_1, Tl_1])
    CRN.write_sbml_file(f'ctRSD_half_and_half.xml') #saving CRN as sbml
    
    with open('temp_CRN_EQNs_half_and_half.txt', 'w') as f:
        f.write(CRN.pretty_print(show_rates = True, show_keys = True))
    
    print('CRN Compiled')
    
    timepoints = np.linspace(0, 12000, 500)
    
    x0 = {assembly_ugate_0.dna:12, assembly_ugate_1.dna:12, assembly_input.dna: 24, 
          'protein_ribozyme_0':1, 'protein_ribozyme_1':1}
    
    data = []
    
    for i in range(100):
        R =  CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, 
                                                 stochastic=True)
        data.append(R['protein_CFP_degtagged'].iloc[-1])
    
    y, x = np.histogram(data, bins=200)
    x_mid = 0.5 * (x[1:] + x[:-1])
    
    plt.fill_between(x_mid, offset, y + offset, color=color, alpha=0.3)
    plt.plot(x_mid, y + offset, color=color, lw=1.2, label = 'half and half')
    plt.plot([1, 2e4], [offset, offset], color='k')
    
def rna_se_biocrnpyler_half_and_half_alt(ktl_0, ktl_1, offset, color):
    
    # Parameters (found in ctRSD_simulator_210.py)
    params = {"ktx":0.013, "kdil":0.00075, "kcat":0.25/60, 
              ParameterKey(mechanism=None, part_id='rna_ugate_0', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_gate_0', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_input_gate_0', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_ugate_1', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_gate_1', name='kdil'):0.000375,
              ParameterKey(mechanism=None, part_id='rna_input_gate_1', name='kdil'):0.000375
              } # kdil not from source, value? 
    #should there be a different kdil for dsRNA
    
    krsd = 1e3/1e9
    krev = 270/1e9
    ktl_0 = ktl_0 # not from source, use what value?
    ktl_1 = ktl_1
    
    # Species
    input_rna = Species('input_rna', material_type='rna')
    output_rna = Species('output_rna', material_type='rna')
    
    ugate_0 = Species('ugate_0', material_type='rna')
    ugate_1 = Species('ugate_1', material_type='rna')
    gate_0 = Species('gate_0', material_type='rna')
    gate_1 = Species('gate_1', material_type='rna')
    input_gate_0 = Species('input_gate_0', material_type='rna')
    input_gate_1 = Species('input_gate_1', material_type='rna')
    
    CFP = Species('CFP', material_type='protein', attributes=['degtagged'])
    
    ribozyme_0 = Enzyme('ribozyme_0', substrates = [ugate_0], products = [gate_0]) #not a protein but defaults to protein material_type
    ribozyme_1 = Enzyme('ribozyme_1', substrates = [ugate_1], products = [gate_1])
    
    # DNAassemblies and Mixture 
    assembly_ugate_0 = DNAassembly(name='ugate_0', promoter='strong', rbs=None,
                                       transcript=ugate_0)
    assembly_ugate_1 = DNAassembly(name='ugate_1', promoter='strong', rbs=None,
                                       transcript=ugate_1)
    assembly_input = DNAassembly(name='input', promoter='strong', rbs=None, 
                                 transcript=input_rna)
    
    dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)
    
    global_mechanisms = {"dilution":dilution_mechanism}
    
    M0 = SimpleTxTlExtract(name="txtl", parameters = params, species=[CFP, input_gate_0, output_rna], 
                          global_mechanisms=global_mechanisms,
                          components=[assembly_ugate_0, assembly_input, 
                                      ribozyme_0])
    
    # Strand Displacement and Translation
    SD_0 = Reaction.from_massaction([input_rna, gate_0], [input_gate_0, output_rna], k_forward=krsd, 
                                  k_reverse=krev)
    Tl_0 = Reaction.from_massaction([input_gate_0], [CFP, input_gate_0], k_forward=ktl_0)
    
    # File Creation and Simulation
    CRN0 = M0.compile_crn()
    CRN0.add_reactions([SD_0, Tl_0])
    CRN0.write_sbml_file(f'ctRSD_half_and_half_alt.xml') #saving CRN as sbml
    
    with open('temp_CRN_EQNs_half_and_half_alt.txt', 'w') as f:
        f.write(CRN0.pretty_print(show_rates = True, show_keys = True))
    
    print('CRN Compiled')
    
    timepoints = np.linspace(0, 12000, 500)
    
    x0 = {assembly_ugate_0.dna:24, assembly_input.dna: 24, 'protein_ribozyme_0':1}
    
    data = []
    
    for i in range(50):
        R =  CRN0.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, 
                                                 stochastic=True)
        data.append(R['protein_CFP_degtagged'].iloc[-1])
        
    M1 = SimpleTxTlExtract(name="txtl", parameters = params, species=[CFP, input_gate_1, output_rna], 
                          global_mechanisms=global_mechanisms,
                          components=[assembly_ugate_1, assembly_input, 
                                      ribozyme_1])
    
    # Strand Displacement and Translation
    SD_1 = Reaction.from_massaction([input_rna, gate_1], [input_gate_1, output_rna], k_forward=krsd, 
                                  k_reverse=krev)
    Tl_1 = Reaction.from_massaction([input_gate_1], [CFP, input_gate_1], k_forward=ktl_1)
    
    # File Creation and Simulation
    CRN1 = M1.compile_crn()
    CRN1.add_reactions([SD_1, Tl_1])
    # CRN1.write_sbml_file(f'ctRSD_half_and_half_alt.xml') #saving CRN as sbml
    
    # with open('temp_CRN_EQNs_half_and_half_alt.txt', 'w') as f:
    #     f.write(CRN.pretty_print(show_rates = True, show_keys = True))
    
    print('CRN Compiled')
    
    x0 = {assembly_ugate_1.dna:24, assembly_input.dna: 24, 'protein_ribozyme_1':1}
    
    for i in range(50):
        R =  CRN1.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, 
                                                 stochastic=True)
        data.append(R['protein_CFP_degtagged'].iloc[-1])
        
    
    y, x = np.histogram(data, bins=200)
    x_mid = 0.5 * (x[1:] + x[:-1])
    
    plt.fill_between(x_mid, offset, y + offset, color=color, alpha=0.3)
    plt.plot(x_mid, y + offset, color=color, lw=1.2, label = 'half and half alt')
    plt.plot([1, 2e4], [offset, offset], color='k')

#parameter list for function
ktl_lst = [0.025, 0.0125, 0.005, 0.0025, 0.001, 0.0005]
color_lst = ['#123524', '#1D5638', '#277A4C', '#32A160', '#4DCB87', '#7FEFB2']
offset_lst = np.arange(len(ktl_lst)) * 5

#Initialize Figure
plt.figure()

rna_se_biocrnpyler_half_and_half(ktl_lst[0], ktl_lst[1], 0, '#B24316')

rna_se_biocrnpyler_half_and_half_alt(ktl_lst[0], ktl_lst[1], offset_lst[1], '#B20000')

for j in range(2):
    rna_se_biocrnpyler(ktl_lst[j], offset_lst[j+2], color_lst[j])
    
plt.title('CFP Stochastic Distribution')
plt.xscale('log')
plt.yticks([])
plt.legend()
plt.xlabel('CFP Signal')
plt.ylabel('Cell Count')

plt.show()
