import selfies as sf

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError, ChemMCPToolProcessError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class Smiles2Selfies(BaseTool):
    __version__ = "0.1.0"
    name = "Smiles2Selfies"
    func_name = 'convert_smiles_to_selfies'
    description = "Convert SMILES to SELFIES string."
    implementation_description = "Uses the SELFIES library to convert a SMILES string to its SELFIES representation."
    categories = ["Molecule"]
    tags = ["Name Conversion", "SMILES", "SELFIES"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('selfies', 'str', 'SELFIES string of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'selfies': '[C][C][O]'}},
    ]
    oss_dependencies = [
        ("selfies", "https://github.com/aspuru-guzik-group/selfies", "Apache License 2.0")
    ]
    services_and_software = []

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")
        try:
            selfies = sf.encoder(smiles)
        except KeyboardInterrupt:
            raise
        except:
            raise ChemMCPToolProcessError("Cannot convert the SMILES into SELFIES, possibly because it is not a valid SMILES string.")
        return selfies


if __name__ == "__main__":
    run_mcp_server()
