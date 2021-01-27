import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm

Water_res_name = 'H2O'
Fake_water_res_name = 'h2o'

FF_file_water = 'files/ffxml/pore-spce.xml'
FF_file_fake_water = 'files/ffxml/FF_Fake_SPCE.xml'

Water_mol2_file = 'files/tip3p.mol2'
Fake_water_mol2_file = 'files/fake_tip3p.mol2'

water = mb.load('O', smiles=True)
water.name = Water_res_name
water.energy_minimize(forcefield = FF_file_water , steps=10**9)

fake_water = mb.load('O', smiles=True)
fake_water.name = Fake_water_res_name
fake_water.energy_minimize(forcefield = FF_file_fake_water , steps=10**9)


water_liq_density_kg_per_m_cubed = 950
number_waters = 1000

Na = 6.022 *10**23
Mw_water = 18.01528
cubic_box_length_m = (number_waters / Na * Mw_water * 1000 / water_liq_density_kg_per_m_cubed * 10**24)**(1/3)/10



FF_Dict = {water.name: FF_file_water }

residues_List = [water.name]

Fix_bonds_angles_residues = [ water.name]



print('Running: filling liquid box')
water_box_liq = mb.fill_box(compound=[water],
                            n_compounds = [number_waters],
                            box=[cubic_box_length_m, cubic_box_length_m, cubic_box_length_m] )
print('Completed: filling liquid box')



print('Running: GOMC FF file, and the psf and pdb files for 600A box')
mf_charmm.charmm_psf_psb_FF(water_box_liq,
                            'Box_0_liq_water_box',
                            structure_1 = None,
                            filename_1 = None,
                            FF_filename ="GOMC_water_FF" ,
                            forcefield_selection = FF_Dict,
                            residues= residues_List ,
                            bead_to_atom_name_dict = None,
                            fix_residue = None,
                            fix_res_bonds_angles = Fix_bonds_angles_residues,
                            reorder_res_in_pdb_psf = False
                            )
print('Completed: GOMC FF file, and the psf and pdb files for 600A box')

