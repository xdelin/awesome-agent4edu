import logging

import synspace
from rdkit import Chem

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.smiles import is_smiles


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MoleculeModifier(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeModifier"
    func_name = 'modify_molecule'
    description = "Proposes small, chemically accessible modifications to a compound using well-established medicinal chemistry reactions and purchasable building blocks. Important: This tool cannot be used more than 3 times in a row or more than 10 times in the entire trajectory."
    implementation_description = "Uses the [synspace](https://github.com/whitead/synspace) package to propose small, chemically accessible modifications to a compound."
    categories = ["Molecule"]
    tags = ["Molecule Modification"]
    required_envs = []
    text_input_sig = [("smiles", "str", "N/A", "The SMILES of the molecule to modify.")]
    code_input_sig = [("smiles", "str", "N/A", "The SMILES of the molecule to modify.")]
    output_sig = [("smiles", "str", "The SMILES of the modified molecule.")]
    examples = [
        {'text_input': {'smiles': 'CCO'}, 'code_input': {'smiles': 'CCO'}, 'output': {'smiles': 'The molecule leads to the following modified molecule(s):\n- C[CH]O: Obtained by a retro Mitsunobu phenole reaction.\n- C[CH]OCCCC1(O)CC1: Obtained by a retro Mitsunobu phenole;Williamson ether reaction.\n'}},
    ]
    oss_dependencies = [
        ("synspace", "https://github.com/whitead/synspace", "MIT")
    ]
    services_and_software = []

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError(f"Invalid SMILES: {smiles}.")
        
        mols, props = synspace.chemical_space(smiles)
        assert len(mols) > 1, "No modifications found."

        text = "The molecule leads to the following modified molecule(s):\n"
        for mol, prop in zip(mols[1:], props[1:]):
            mol_smiles = Chem.MolToSmiles(mol)
            rxn_name = prop['rxn-name']
            rxn_name = rxn_name.replace('-', ' ')
            text += f"- {mol_smiles}: Obtained by a {rxn_name} reaction.\n"

        return text


if __name__ == "__main__":
    run_mcp_server()
