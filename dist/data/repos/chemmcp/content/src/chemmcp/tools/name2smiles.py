import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPSearchFailError
from ..tool_utils.names import pubchem_name2smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Name2Smiles(BaseTool):
    __version__ = "0.1.0"
    name = "Name2Smiles"
    func_name = 'convert_chemical_name_to_smiles'
    description = "Convert chemical name to SMILES string."
    implementation_description = "Uses PubChem's API to convert a chemical name to its corresponding SMILES representation. The conversion is performed by searching PubChem's database for the molecule and retrieving its SMILES string."
    categories = ["Molecule"]
    tags = ["Name Conversion", "SMILES", "Molecular Names", "PubChem", "APIs"]
    required_envs = []
    code_input_sig = [('name', 'str', 'N/A', 'Chemical name of the molecule')]
    text_input_sig = [('name', 'str', 'N/A', 'Chemical name of the molecule')]
    output_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    examples = [
        {'code_input': {'name': 'aspirin'}, 'text_input': {'name': 'aspirin'}, 'output': {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'}},
    ]
    oss_dependencies = [
        ("PubChemPy", "https://github.com/mcs07/PubChemPy", "MIT")
    ]
    services_and_software = []

    def _run_base(self, name: str) -> str:
        try:
            smi = pubchem_name2smiles(name)
            logger.debug("Looking up PubChem succeeded.")
        except ChemMCPSearchFailError as e:
            logger.debug("Looking up PubChem failed.")
            raise e
        
        return smi


if __name__ == "__main__":
    run_mcp_server()
