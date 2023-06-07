# Reorder chains and renumber residues of AlphaFold2 models
# of the tetramer-of-heterodimers of NylC homologs.
#
# JMP 3/23/2023
# 
# Usage:
# python post-process_af2_model.py 

from pymol import cmd
import argparse
import sys
import json
import numpy

def load_pdb_files():
  cmd.load(ref)
  cmd.load(model)
  objs = cmd.get_object_list("(all)")

  cmd.set_name(objs[0], ref_name)
  cmd.set_name(objs[1], model_name)

  cmd.hide("everything", "nylc and resn HOH+PO4+CL+GOL")
  cmd.set_view ((
    -0.033802517,    0.158777967,    0.986735702,
    -0.999350965,   -0.017669253,   -0.031391911,
     0.012449547,   -0.987155139,    0.159273610,
    -0.000052944,    0.000145555, -306.412078857,
   -21.511169434,  -33.441307068,   42.405929565,
   253.244506836,  359.581939697,  -20.000000000 ))

  cmd.util.cbc()

def get_model_chain_ids():
  """Get the chain IDs in the AF2 model."""
  model_chains = cmd.get_chains(model_name)

  # The alpha chains should always be the first four and the beta chains should always be the last four.

  alpha_chains = model_chains[:4]
  beta_chains = model_chains[-4:]

  return alpha_chains, beta_chains

def align_model_to_ref():
  """Align the first chain in the AF2 structure to chain A of NylC"""
  ref_selection = ref_name + ' and chain A'
  model_selection = model_name + ' and chain ' + alpha_chains[0]
  cmd.align(model_selection, ref_selection)
  #cmd.super(model_selection, ref_selection)

def assign_alpha_chains_temp_names():
  """Assign alpha chains temporary chain IDs (W, X, Y, and Z) to avoid accidental overwrites"""

  for i, chain in enumerate(alpha_chains):
    alter_from = '(' + model_name + ' and chain ' + chain +')'
    alter_to = 'chain="' + tmp_chains[i] +'"'
    cmd.alter(alter_from, alter_to)

  cmd.sort(model_name)

def update_alpha_chain_order():
  """Reorder the alpha chains to match the NylC X-ray template"""
  ref_chains = cmd.get_chains(ref_name)
  model_chains = cmd.get_chains(model_name)
  
  model_alpha_chains = model_chains[-4:]
  
  new_alpha_chains = []
  for i in ref_chains:
    rmsd_list = []
    isel = ref_name + ' and chain ' + i + ' and resi 5-260' 
    isel2 = isel + ' & aln'
    
    for j in model_alpha_chains:
      jsel = model_name + ' and chain ' + j
      jsel2 = jsel + ' & aln'
      cmd.align(isel, jsel, cycles=0, transform=0, object="aln")

      # check alignment
      #raw_aln = cmd.get_raw_alignment('aln')
      #idx2resi = {}
      #cmd.iterate('aln', 'idx2resi[model, index] = resi', space={'idx2resi': idx2resi})
      #for idx1, idx2 in raw_aln:
      #  print('%s -> %s' % (idx2resi[idx1], idx2resi[idx2]))
      #

      rms_cur = cmd.rms_cur(jsel2, isel2, matchmaker=-1)
      rmsd_list.append(rms_cur)
      cmd.delete("aln")

    # Get the index of the lowest RMS value in each list. That's the chain that should be renamed to
    # match the chain we are currently considering.

    min_rms = min(rmsd_list)
    min_rms_index = rmsd_list.index(min_rms)
    new_alpha_chains.append(model_alpha_chains[min_rms_index])

  return new_alpha_chains

def rename_chains(new_chains):
  for i, chain in enumerate(new_chains):
    alter_from = '(' + model_name + ' and chain ' + chain +')'
    alter_to = 'chain="' + alpha_chains[i] +'"'
    cmd.alter(alter_from, alter_to)

  cmd.sort(model_name)
  cmd.util.cbc()

def renumber_beta_chains():

  alter_from = '(' + model_name + ' and chain E+F+G+H)' 
  offset = model_beta_resi_start-1
  alter_to = 'resi=str(int(resi)+'+str(offset)+')'
  cmd.alter(alter_from, alter_to) 
  cmd.sort()

def update_beta_chain_order():
  """Reorder the beta chains to match the NylC X-ray template"""
  ref_chains = cmd.get_chains(ref_name)
  model_chains = cmd.get_chains(model_name)
  
  model_beta_chains = model_chains[-4:]
  
  new_beta_chains = []
  for i in ref_chains:
    rmsd_list = []
    isel = ref_name + ' and chain ' + i + ' and resi 267-340' 
    isel2 = isel + ' & aln'
    
    for j in model_beta_chains:
      jsel = model_name + ' and chain ' + j
      jsel2 = jsel + ' & aln'
      cmd.align(isel, jsel, cycles=0, transform=0, object="aln")
      rms_cur = cmd.rms_cur(jsel2, isel2, matchmaker=-1)
      rmsd_list.append(rms_cur)
      cmd.delete("aln")

    min_rms = min(rmsd_list)
    min_rms_index = rmsd_list.index(min_rms)
    new_beta_chains.append(model_beta_chains[min_rms_index])

  return new_beta_chains

def write_pdb():
  #pdbname = model_name + '_fixed.pdb'
  pdbname = model_name + '.pdb'
  # Preserve the coordinates from the alignment to the template when writing the new PDB file.
  m = cmd.get_view(1)
  ttt = [m[0], m[1], m[2], 0.0,
         m[3], m[4], m[5], 0.0,
         m[6], m[7], m[8], 0.0,
         0.0,   0.0,  0.0, 1.0]
  cmd.transform_object(object,ttt,transpose=1)
  cmd.save(pdbname, model_name)

# MAIN

# These lines should not be changed
ref="../5y0m_tetramer.pdb"
#ref="5y0m_tetramer.pdb"
ref_name = "nylc"
#
# Get command line arguments.
parser = argparse.ArgumentParser()
parser.add_argument("--model", help='AlphaFold2 rank 1 model (PDB)', type=str)
parser.add_argument("--beta_start", help='residue number of beta chain N-terminus', type=int)
args = parser.parse_args()

model=args.model
model_name=args.model.split('_4x2')[0]
model_beta_resi_start=args.beta_start # Offset beta chain residue numbers by this number minus one.
tmp_chains = [ 'W', 'X', 'Y', 'Z']

load_pdb_files()
alpha_chains, beta_chains = get_model_chain_ids()
align_model_to_ref()
assign_alpha_chains_temp_names()
new_alpha_chains = update_alpha_chain_order()
rename_chains(new_alpha_chains)
renumber_beta_chains()
new_beta_chains = update_beta_chain_order()
rename_chains(new_beta_chains)
write_pdb()
