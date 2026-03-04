import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPToolProcessError
from ..tool_utils.smiles import is_smiles
from ..tool_utils.names import smiles2formula
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Smiles2Formula(BaseTool):
    __version__ = "0.1.0"
    name = "Smiles2Formula"
    func_name = 'convert_smiles_to_formula'
    description = "Convert SMILES to molecular formula."
    implementation_description = "Uses RDKit to parse the SMILES string and generate the corresponding molecular formula. The formula is calculated by counting atoms and their occurrences in the molecule."
    oss_dependencies = []
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]
    categories = ["Molecule"]
    tags = ["Name Conversion", "Molecular Formulas", "SMILES", "RDKit"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('formula', 'str', 'Molecular formula of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'formula': 'C2H6O'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")
        
        try:
            formula = smiles2formula(smiles)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMCPToolProcessError("Failed to process the SMILES string.") from e
        
        return formula


if __name__ == "__main__":
    run_mcp_server()
