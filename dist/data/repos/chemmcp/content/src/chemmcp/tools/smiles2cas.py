import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPSearchFailError
from ..tool_utils.smiles import is_smiles
from ..tool_utils.names import pubchem_smiles2cas
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Smiles2Cas(BaseTool):
    __version__ = "0.1.0"
    name = "Smiles2Cas"
    func_name = 'convert_smiles_to_cas'
    description = "Convert SMILES to CAS number based on PubChem."
    implementation_description = "Searches PubChem for the molecule and returns the CAS number."
    categories = ["Molecule"]
    tags = ["Name Conversion", "SMILES", "CAS", "PubChem", "APIs"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('cas', 'str', 'CAS number of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'cas': '64-17-5'}},
    ]
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT")
    ]
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")
        
        try:
            cas = pubchem_smiles2cas(smiles)
            logger.debug("Looking up PubChem succeeded.")
        except ChemMCPSearchFailError as e:
            logger.debug("Looking up PubChem failed.")
            raise e
        
        return cas


if __name__ == "__main__":
    run_mcp_server()
