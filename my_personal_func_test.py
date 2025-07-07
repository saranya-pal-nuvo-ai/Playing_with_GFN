from mols.mol_mdp_ext import BlockMoleculeData
from rdkit import Chem
from mols.utils.chem import mol_from_frag  # assuming it's the function used
from mols.utils.molMDP import BlockMoleculeData
from mols.mol_mdp_ext import BlockMoleculeDataExtended
from mols.mol_mdp_ext import MolMDPExtended

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



def test_Mol_MDP(json_file):
    N = 3

    mdp = MolMDP(json_file)

    # print(mdp.block_smi[:N])
    # print(mdp.block_rs[:N])
    # print(mdp.block_nrs[:N])
    # print(mdp.block_natm[:N])

    
    # mdp.add_block(
    #     block_idx=,
    #     stem_idx=,
    #     atmidx=,
        
    # )



def a_complete_dry_run(json_file):
    # 1) load the fragment library
    mdp = MolMDPExtended(json_file)
    mdp.post_init(device=None, repr_type="block_graph")
    mdp.build_translation_table()
    mdp.reset()   # start from the empty molecule


    # 2) find the library‐indices for “C” and “O”
    #    (they may not be the first entries in mdp.block_smi; let's search)
    for i, smi in enumerate(mdp.block_smi):
        if smi in ("C","O"):
            print(i, smi, mdp.block_rs[i])
    # Suppose this prints:
    #   42  C  [0]        ← a single‐atom carbon block with one R‐site
    #   17  O  [0]        ← a single‐atom oxygen block

    # 3) add the methyl (C) fragment as the first block:
    c_idx = mdp.block_smi.index("C")   # e.g. 42
    mdp.add_block(block_idx=c_idx)     # stem_idx=None for the very first block

    print("After adding C:")
    print("  blockidxs:",   mdp.molecule.blockidxs)       # [42]
    print("  stems:",       mdp.molecule.stems)           # e.g. [[0,0]]
    print("  smiles:",      mdp.molecule.smiles)          # "C"

    # 4) now attach the hydroxyl (O) to the only open stem:
    o_idx = mdp.block_smi.index("O")   # e.g. 17
    stem_to_use = 0                     # the first (and only) stem in mdp.molecule.stems
    mdp.add_block(block_idx=o_idx, stem_idx=stem_to_use)

    print("After adding O onto that stem:")
    print("  blockidxs:",   mdp.molecule.blockidxs)      # [42, 17]
    print("  jbonds:",      mdp.molecule.jbonds)        # [[0,1, 0,0]]  C(0)–O(0)
    print("  stems:",       mdp.molecule.stems)         # [] (no more open sites)
    print("  smiles:",      mdp.molecule.smiles)        # "CO" (methanol)
    


def main():
    test_functions()


if __name__ == '__main__':
    # main()

    # test_mdp()
    json_file = '{"block_name":{"0":"c1ccccc1_0","1":"CO_0","2":"CO_1","3":"C=O_0","4":"C1CCNCC1_0","5":"C1CCNCC1_1","6":"C1CCNCC1_5","7":"O_0","8":"C_0","9":"c1ccncc1_0","10":"c1ccncc1_1","11":"c1ccncc1_3","12":"CNC=O_0","13":"CNC=O_2","14":"CC_0","15":"O=CNO_2","16":"O=CO_2","17":"CCC_2","18":"c1ccc2ccccc2c1_0","19":"c1ccc2ccccc2c1_4","20":"NC=O_0","21":"NC=O_2","22":"C1CCCC1_0","23":"CN_0","24":"CN_1","25":"F_0","26":"N_0","27":"CNC_2","28":"O=PO_2","29":"Cl_0","30":"C1CCOC1_0","31":"C1CCOC1_2","32":"c1ncc2nc[nH]c2n1_2","33":"c1ncc2nc[nH]c2n1_6","34":"c1cc[nH]c1_0","35":"c1cc[nH]c1_2","36":"c1ccc2[nH]ccc2c1_4","37":"C1COCCN1_4","38":"C1CCCCC1_0","39":"N=CN_2","40":"c1cscn1_0","41":"c1cscn1_2","42":"O=S=O_2","43":"C1CCNC1_0","44":"C1CCNC1_2","45":"C1CCNC1_4","46":"c1cn[nH]c1_0","47":"c1cn[nH]c1_1","48":"c1cn[nH]c1_2","49":"c1cn[nH]c1_4","50":"O=[NH+][O-]_2","51":"FC(F)F_3","52":"c1ccsc1_2","53":"O=[PH](O)O_3","54":"S_0","55":"CC(C)C_3","56":"I_0","57":"O=c1[nH]cnc2[nH]cnc12_1","58":"O=c1[nH]cnc2[nH]cnc12_6","59":"CNC(C)=O_0","60":"c1cncnc1_0","61":"c1cncnc1_1","62":"c1cncnc1_3","63":"Br_0","64":"C1CNCCN1_4","65":"c1cc[nH+]cc1_1","66":"c1cc[nH+]cc1_5","67":"C1CC1_0","68":"C=C_0","69":"C1CCOCC1_0","70":"C1CCOCC1_1","71":"C1CCOCC1_3","72":"CC(C)O_3","73":"N[SH](=O)=O_3","74":"CCO_2","75":"CC(N)=O_1","76":"C#N_0","77":"CC=O_2","78":"O=c1cc[nH]c(=O)[nH]1_3","79":"O=c1cc[nH]c(=O)[nH]1_4","80":"C1=CNC=CC1_0","81":"C1=CNC=CC1_5","82":"C=CC_0","83":"C=CC_2","84":"CS_0","85":"O=[SH](=O)O_3","86":"C[NH3+]_0","87":"C=C(C)C_0","88":"C=N_0","89":"O=C[O-]_2","90":"O=c1nccc[nH]1_2","91":"O=c1nccc[nH]1_5","92":"O=[SH](=O)[O-]_3","93":"C1=CCCCC1_0","94":"C1=CCCCC1_2","95":"C1=CCCCC1_4","96":"O=c1[nH]cnc2c1NCCN2_1","97":"O=c1[nH]cnc2c1NCCN2_3","98":"O=P[O-]_2","99":"O=[PH]([O-])O_3","100":"O=c1nc2[nH]c3ccccc3nc-2c(=O)[nH]1_2","101":"O=c1nc2[nH]c3ccccc3nc-2c(=O)[nH]1_3","102":"O=c1nc2[nH]c3ccccc3nc-2c(=O)[nH]1_9","103":"C[SH2+]_1","104":"C1COCC[NH2+]1_4"},"block_smi":{"0":"c1ccccc1","1":"CO","2":"CO","3":"C=O","4":"C1CCNCC1","5":"C1CCNCC1","6":"C1CCNCC1","7":"O","8":"C","9":"c1ccncc1","10":"c1ccncc1","11":"c1ccncc1","12":"CNC=O","13":"CNC=O","14":"CC","15":"O=CNO","16":"O=CO","17":"CCC","18":"c1ccc2ccccc2c1","19":"c1ccc2ccccc2c1","20":"NC=O","21":"NC=O","22":"C1CCCC1","23":"CN","24":"CN","25":"F","26":"N","27":"CNC","28":"O=PO","29":"Cl","30":"C1CCOC1","31":"C1CCOC1","32":"c1ncc2nc[nH]c2n1","33":"c1ncc2nc[nH]c2n1","34":"c1cc[nH]c1","35":"c1cc[nH]c1","36":"c1ccc2[nH]ccc2c1","37":"C1COCCN1","38":"C1CCCCC1","39":"N=CN","40":"c1cscn1","41":"c1cscn1","42":"O=S=O","43":"C1CCNC1","44":"C1CCNC1","45":"C1CCNC1","46":"c1cn[nH]c1","47":"c1cn[nH]c1","48":"c1cn[nH]c1","49":"c1cn[nH]c1","50":"O=[NH+][O-]","51":"FC(F)F","52":"c1ccsc1","53":"O=[PH](O)O","54":"S","55":"CC(C)C","56":"I","57":"O=c1[nH]cnc2[nH]cnc12","58":"O=c1[nH]cnc2[nH]cnc12","59":"CNC(C)=O","60":"c1cncnc1","61":"c1cncnc1","62":"c1cncnc1","63":"Br","64":"C1CNCCN1","65":"c1cc[nH+]cc1","66":"c1cc[nH+]cc1","67":"C1CC1","68":"C=C","69":"C1CCOCC1","70":"C1CCOCC1","71":"C1CCOCC1","72":"CC(C)O","73":"N[SH](=O)=O","74":"CCO","75":"CC(N)=O","76":"C#N","77":"CC=O","78":"O=c1cc[nH]c(=O)[nH]1","79":"O=c1cc[nH]c(=O)[nH]1","80":"C1=CNC=CC1","81":"C1=CNC=CC1","82":"C=CC","83":"C=CC","84":"CS","85":"O=[SH](=O)O","86":"C[NH3+]","87":"C=C(C)C","88":"C=N","89":"O=C[O-]","90":"O=c1nccc[nH]1","91":"O=c1nccc[nH]1","92":"O=[SH](=O)[O-]","93":"C1=CCCCC1","94":"C1=CCCCC1","95":"C1=CCCCC1","96":"O=c1[nH]cnc2c1NCCN2","97":"O=c1[nH]cnc2c1NCCN2","98":"O=P[O-]","99":"O=[PH]([O-])O","100":"O=c1nc2[nH]c3ccccc3nc-2c(=O)[nH]1","101":"O=c1nc2[nH]c3ccccc3nc-2c(=O)[nH]1","102":"O=c1nc2[nH]c3ccccc3nc-2c(=O)[nH]1","103":"C[SH2+]","104":"C1COCC[NH2+]1"},"block_r":{"0":[0,1,2,3,4,5],"1":[0,1],"2":[1,0,0,0],"3":[0],"4":[0,3],"5":[1,0,3],"6":[3,0],"7":[0],"8":[0],"9":[0,1,2,4,5],"10":[1,0,2,4,5],"11":[2,0,1,4,5],"12":[0,2],"13":[2,0,0],"14":[0,1],"15":[1],"16":[1],"17":[1],"18":[0],"19":[2],"20":[0,1],"21":[1,0],"22":[0],"23":[0,1,1],"24":[1,0,0],"25":[0],"26":[0],"27":[1],"28":[1],"29":[0],"30":[0,4,1,2,4],"31":[2,4,0,1,4],"32":[2,6],"33":[6,2],"34":[0],"35":[2],"36":[6],"37":[5],"38":[0,1,2,3,4,5],"39":[1],"40":[0,3],"41":[3,0],"42":[1],"43":[0,3,4],"44":[2,0,3,4],"45":[3,0,4],"46":[0,1,3,4],"47":[1,0,3,4],"48":[4,0,1,3],"49":[3,0,1,4],"50":[1],"51":[1],"52":[2,4],"53":[1],"54":[0],"55":[1],"56":[0],"57":[3,6],"58":[6,3],"59":[0],"60":[0,1,3,5],"61":[1,0,3,5],"62":[3,0,1,5],"63":[0],"64":[2,5],"65":[1,3],"66":[3,1],"67":[0],"68":[0,1],"69":[0,1,2,4,5],"70":[1,0,2,4,5],"71":[2,0,1,4,5],"72":[1],"73":[1],"74":[1],"75":[2],"76":[0],"77":[1],"78":[2,4],"79":[4,2],"80":[0,2],"81":[2,0],"82":[0,1],"83":[1,0],"84":[0],"85":[1],"86":[0],"87":[0],"88":[0],"89":[1],"90":[3,6],"91":[6,3],"92":[1],"93":[0],"94":[2,0],"95":[3,0],"96":[3,8],"97":[8,3],"98":[1],"99":[1],"100":[8,4,7],"101":[7,4,8],"102":[4,7,8],"103":[1],"104":[5]}}'
    # test_Mol_MDP(json_file)
    a_complete_dry_run(json_file)