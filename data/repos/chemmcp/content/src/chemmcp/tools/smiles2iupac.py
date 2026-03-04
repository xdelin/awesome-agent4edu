import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPSearchFailError
from ..tool_utils.smiles import is_smiles
from ..tool_utils.names import pubchem_smiles2iupac
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Smiles2Iupac(BaseTool):
    __version__ = "0.1.0"
    name = "Smiles2Iupac"
    func_name = 'convert_smiles_to_iupac'
    description = "Convert SMILES to IUPAC name."
    implementation_description = "Uses PubChem's API to convert a SMILES string to its corresponding IUPAC name. The conversion is performed by searching PubChem's database for the molecule and retrieving its IUPAC name."
    categories = ["Molecule"]
    tags = ["Name Conversion", "SMILES", "IUPAC", "PubChem", "APIs"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('iupac', 'str', 'IUPAC name of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'iupac': 'ethanol'}},
    ]
    oss_dependencies = [
        ("PubChemPy", "https://github.com/mcs07/PubChemPy", "MIT")
    ]
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")
        
        try:
            name = pubchem_smiles2iupac(smiles)
            logger.debug("Looking up PubChem succeeded.")
        except KeyboardInterrupt:
            raise
        except ChemMCPSearchFailError as e:
            logger.debug("Looking up PubChem failed.")
            raise e
        
        # If PubChem fails, try STOUT
        # TODO: The STOUT package is not available anymore. Should be fixed later.
        
        return name


if __name__ == "__main__":
    run_mcp_server()
