# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 19:27:04 2025

@author: zacha
"""

from biocrnpyler import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters
params = {"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.05, "kdeg":0.001, 
          "kdil":0.0075, "kcat":10}
complex_params = {"ku":20, "kb":100}

# Species
input_rna = Species('input_rna', material_type='rna')
output_rna = Species('output_rna', material_type='rna')
ugate = Species('ugate', material_type='rna')
gate = Species('gate', material_type='rna')

input_gate = ChemicalComplex([input_rna, gate], 'input_gate', material_type='rna', 
                             parameters=complex_params)

CFP = Species('CFP', material_type='protein')

ribozyme = Enzyme('ribozyme', substrates = [ugate], products = [gate]) #not a protein but defaults to protein material_type

# DNAassemblies and Mixture 
assembly_ugate = DNAassembly(name='ugate', promoter='strong', rbs=None,
                                   transcript=ugate)
assembly_input = DNAassembly(name='input', promoter='strong', rbs=None, 
                             transcript=input_rna)

M = SimpleTxTlExtract(name="txtl", parameters = params,
                      components=[assembly_ugate, assembly_input, input_gate, 
                                  ribozyme])

# File Creation and Simulation
CRN = M.compile_crn()
CRN.write_sbml_file('rna-strand-exchange.xml') #saving CRN as sbml

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

print('CRN Compiled')

timepoints = np.linspace(0, 1000, 200)

x0 = {assembly_ugate.dna:1, assembly_input.dna: 1, 'protein_ribozyme':1}

R = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, safe=True)

# Plotting
plot_lst = ['rna_ugate', 'rna_gate', 'rna_input_rna', 'rna_input_gate_']

plt.figure()

for i in range(len(plot_lst)):
    plt.plot(R['time'], R[plot_lst[i]], label = plot_lst[i])

plt.legend()
plt.title('')
plt.xlabel('')
plt.ylabel('')

plt.show()
