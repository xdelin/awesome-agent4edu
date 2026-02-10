from rdkit import Chem

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


# Constants
# List obtained from https://github.com/rdkit/rdkit/blob/master/Data/FunctionalGroups.txt
DICT_FGS = {
    "furan": "o1cccc1",
    "aldehydes": " [CX3H1](=O)[#6]",
    "esters": " [#6][CX3](=O)[OX2H0][#6]",
    "ketones": " [#6][CX3](=O)[#6]",
    "amides": " C(=O)-N",
    "thiol groups": " [SH]",
    "alcohol groups": " [OH]",
    "methylamide": "*-[N;D2]-[C;D3](=O)-[C;D1;H3]",
    "carboxylic acids": "*-C(=O)[O;D1]",
    "carbonyl methylester": "*-C(=O)[O;D2]-[C;D1;H3]",
    "terminal aldehyde": "*-C(=O)-[C;D1]",
    "amide": "*-C(=O)-[N;D1]",
    "carbonyl methyl": "*-C(=O)-[C;D1;H3]",
    "isocyanate": "*-[N;D2]=[C;D2]=[O;D1]",
    "isothiocyanate": "*-[N;D2]=[C;D2]=[S;D1]",
    "nitro": "*-[N;D3](=[O;D1])[O;D1]",
    "nitroso": "*-[N;R0]=[O;D1]",
    "oximes": "*=[N;R0]-[O;D1]",
    "Imines": "*-[N;R0]=[C;D1;H2]",
    "terminal azo": "*-[N;D2]=[N;D2]-[C;D1;H3]",
    "hydrazines": "*-[N;D2]=[N;D1]",
    "diazo": "*-[N;D2]#[N;D1]",
    "cyano": "*-[C;D2]#[N;D1]",
    "primary sulfonamide": "*-[S;D4](=[O;D1])(=[O;D1])-[N;D1]",
    "methyl sulfonamide": "*-[N;D2]-[S;D4](=[O;D1])(=[O;D1])-[C;D1;H3]",
    "sulfonic acid": "*-[S;D4](=O)(=O)-[O;D1]",
    "methyl ester sulfonyl": "*-[S;D4](=O)(=O)-[O;D2]-[C;D1;H3]",
    "methyl sulfonyl": "*-[S;D4](=O)(=O)-[C;D1;H3]",
    "sulfonyl chloride": "*-[S;D4](=O)(=O)-[Cl]",
    "methyl sulfinyl": "*-[S;D3](=O)-[C;D1]",
    "methyl thio": "*-[S;D2]-[C;D1;H3]",
    "thiols": "*-[S;D1]",
    "thio carbonyls": "*=[S;D1]",
    "halogens": "*-[#9,#17,#35,#53]",
    "t-butyl": "*-[C;D4]([C;D1])([C;D1])-[C;D1]",
    "tri fluoromethyl": "*-[C;D4](F)(F)F",
    "acetylenes": "*-[C;D2]#[C;D1;H]",
    "cyclopropyl": "*-[C;D3]1-[C;D2]-[C;D2]1",
    "ethoxy": "*-[O;D2]-[C;D2]-[C;D1;H3]",
    "methoxy": "*-[O;D2]-[C;D1;H3]",
    "side-chain hydroxyls": "*-[O;D1]",
    "ketones": "*=[O;D1]",
    "primary amines": "*-[N;D1]",
    "nitriles": "*#[N;D1]",
}
    

@ChemMCPManager.register_tool
class FunctionalGroups(BaseTool):
    __version__ = "0.1.0"
    name = "FunctionalGroups"
    func_name = 'identify_functional_groups'
    description = "Identify functional groups in a molecule."
    implementation_description = "Uses RDKit's SMARTS patterns to identify functional groups in a molecule. The tool checks for a comprehensive list of common functional groups including alcohols, aldehydes, ketones, carboxylic acids, and many others."
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT")
    ]
    services_and_software = []
    categories = ["Molecule"]
    tags = ["Molecular Information", "RDKit", "SMILES"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('fgs', 'str', 'A description of functional groups in the molecule.')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'fgs': 'This molecule contains alcohol groups, and side-chain hydroxyls.'}},
    ]

    def _run_base(self, smiles: str) -> str:
        def _is_fg_in_mol(mol, fg):
            fgmol = Chem.MolFromSmarts(fg)
            mol = Chem.MolFromSmiles(mol.strip())
            return len(Chem.Mol.GetSubstructMatches(mol, fgmol, uniquify=True)) > 0
        
        try:
            fgs_in_molec = [
                name
                for name, fg in DICT_FGS.items()
                if _is_fg_in_mol(smiles, fg)
            ]
            if len(fgs_in_molec) > 1:
                return f"This molecule contains {', '.join(fgs_in_molec[:-1])}, and {fgs_in_molec[-1]}."
            elif len(fgs_in_molec) == 1:
                return f"This molecule contains {fgs_in_molec[0]}."
            else:
                return "This molecule does not contain any functional groups."
        except:
            raise ChemMCPInputError("Invalid SMILES string.")


if __name__ == "__main__":
    run_mcp_server()

