import selfies as sf

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPToolProcessError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class Selfies2Smiles(BaseTool):
    __version__ = "0.1.0"
    name = "Selfies2Smiles"
    func_name = 'convert_selfies_to_smiles'
    description = "Convert SELFIES to SMILES string."
    implementation_description = "Uses the SELFIES library to convert a SELFIES string back to its SMILES representation."
    categories = ["Molecule"]
    tags = ["SMILES", "SELFIES", "Name Conversion"]
    required_envs = []
    code_input_sig = [('selfies', 'str', 'N/A', 'SELFIES string of the molecule')]
    text_input_sig = [('selfies', 'str', 'N/A', 'SELFIES string of the molecule')]
    output_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    examples = [
        {'code_input': {'selfies': '[C][C][O]'}, 'text_input': {'selfies': '[C][C][O]'}, 'output': {'smiles': 'CCO'}},
    ]
    oss_dependencies = [
        ("selfies", "https://github.com/aspuru-guzik-group/selfies", "Apache License 2.0")
    ]
    services_and_software = []

    def _run_base(self, selfies: str) -> str:
        try:
            smiles = sf.decoder(selfies)
        except KeyboardInterrupt:
            raise
        except:
            raise ChemMCPToolProcessError("Cannot convert the SELFIES into SMILES, possibly because it is not a valid SELFIES string.")
        return smiles


if __name__ == "__main__":
    run_mcp_server()
