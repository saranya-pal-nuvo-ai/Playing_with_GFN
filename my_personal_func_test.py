from mols.mol_mdp_ext import BlockMoleculeData
from rdkit import Chem
from mols.utils.chem import mol_from_frag  # assuming it's the function used
from mols.utils.molMDP import BlockMoleculeData
from mols.mol_mdp_ext import BlockMoleculeDataExtended

from mols.utils.molMDP import BlockMoleculeData, MolMDP


def test_functions():
    frag1 = Chem.MolFromSmiles("C[*]")  # C-R

    # Fragment 2: Ethyl group with attachment point
    frag2 = Chem.MolFromSmiles("*CC")  # R-CH2-CH3

    # Join bond between frag1 (atom 1 = *) and frag2 (atom 0 = *)
    jbonds = [(0, 1, 1, 0)]  # (frag1_idx, atom_idx, frag2_idx, atom_idx)

    moldata = BlockMoleculeDataExtended()
    moldata.blocks = [frag1, frag2]
    moldata.blockidxs = [0, 1]
    moldata.jbonds = jbonds
    moldata.numblocks = 2
    moldata.slices = [0, 1]
    moldata.stems = []

    print("SMILES:", moldata.smiles)



def test_mdp():
    """
    We will try to build CH3-OH using 2 blocks: CH3 and OH
    """
    me_alc = BlockMoleculeData()   # Initially molecule is empty

    # Add methyl
    me_alc.add_block(
            block_idx = 0,
            block = Chem.MolFromSmiles("C"),
            block_r = [0],
            stem_idx = None,
            atmidx = None
        )

    print(me_alc.stems)
    
    # Add hydroxide
    me_alc.add_block(
        block_idx=1,
        block=Chem.MolFromSmiles("O"),
        block_r=[0],
        stem_idx=0,
        atmidx=None
    )


    print('-'*10)
    print(me_alc.slices)
    print(me_alc.jbonds)
    print(me_alc.stems)
    


def main():
    test_functions()


if __name__ == '__main__':
    # main()

    test_mdp()