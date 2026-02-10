import logging

from ..utils.base_tool import BaseTool
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.smiles import is_smiles
from ..utils.errors import ChemMCPInputError


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class MoleculeSmilesCheck(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeSmilesCheck"
    func_name = 'check_molecule_smiles'
    description = "Check the syntactical validity of a molecular SMILES string."
    implementation_description = "Uses the [RDKit](https://www.rdkit.org/) library to check the syntactical validity of a molecular SMILES string. Uses a text template to construct textual output."
    oss_dependencies = []
    services_and_software = []
    categories = ["Molecule"]
    tags = ["SMILES", "RDKit", "Molecular Information"]
    required_envs = []
    text_input_sig = [("smiles", "str", "N/A", "The SMILES string to check.")]
    code_input_sig = [("smiles", "str", "N/A", "The SMILES string to check.")]
    output_sig = [("result", "str", "Description of the validity of the SMILES string.")]
    examples = [
        {
            'text_input': {'smiles': 'CCO'}, 
            'code_input': {'smiles': 'CCO'}, 
            'output': {'result': 'The SMILES string is valid.'}
        },
    ]

    def _run_base(self, smiles: str) -> str:
        if '>' in smiles:
            raise ChemMCPInputError("The input contains \">\", which indicates that it may be a reaction SMARTS instead of a molecular SMILES. Please use the SmartsCheck tool instead.")
        
        if is_smiles(smiles):
            return "The molecular SMILES string is valid."
        else:
            return "The molecular SMILES string is invalid."


if __name__ == "__main__":
    run_mcp_server()
