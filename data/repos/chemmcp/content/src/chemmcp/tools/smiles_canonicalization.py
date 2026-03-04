from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.canonicalization import canonicalize_molecule_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class SmilesCanonicalization(BaseTool):
    __version__ = "0.1.0"
    name = "SmilesCanonicalization"
    func_name = 'canonicalize_smiles'
    description = "Canonicalize a molecular SMILES string."
    implementation_description = "Uses a customized version of RDKit's SMILES canonicalization functionality to convert a SMILES string into its canonical form."
    categories = ["Molecule"]
    tags = ["SMILES", "RDKit", "Molecular Operations"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.'), ('isomeric', 'bool', 'N/A', 'Whether to include isomeric information. Default is True.'), ('kekulization', 'bool', 'N/A', 'Whether to perform kekulization. Default is True.'), ('keep_atom_map', 'bool', 'N/A', 'Whether to keep atom mapping numbers, if any. Default is True.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]  # TODO: Support options.
    output_sig = [('canonical_smiles', 'str', 'Canonicalized SMILES string.')]
    examples = [
        {'code_input': {'smiles': 'C(O)C', 'isomeric': True, 'kekulization': True, 'keep_atom_map': False}, 'text_input': {'smiles': 'C(O)C'}, 'output': {'canonical_smiles': 'CCO'}},
    ]
    oss_dependencies = [
        ("LlaSMol", "https://github.com/OSU-NLP-Group/LLM4Chem", "MIT")
    ]
    services_and_software = []

    def _run_base(self, smiles: str, isomeric: bool = True, kekulization: bool = True, keep_atom_map: bool = True) -> str:
        smiles = canonicalize_molecule_smiles(smiles, isomeric=isomeric, kekulization=kekulization, keep_atom_map=keep_atom_map)
        if smiles is None:
            raise ChemMCPInputError("Invalid SMILES string.")
        return smiles
    
    def _run_text(self, smiles: str) -> str:
        return self._run_base(smiles)


if __name__ == "__main__":
    run_mcp_server()

