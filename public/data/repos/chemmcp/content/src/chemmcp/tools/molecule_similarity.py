from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..tool_utils.smiles import tanimoto, is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class MoleculeSimilarity(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeSimilarity"
    func_name = 'cal_molecule_similarity'
    description = "Get the Tanimoto similarity of two molecules. It Can also be used to check if two molecules are identical."
    implementation_description = "Uses the [RDKit](https://www.rdkit.org/) library to calculate the Tanimoto similarity of two molecules. Uses a text template to construct textual output."
    categories = ["Molecule"]
    tags = ["Molecular Information", "RDKit", "SMILES", "Molecular Operations"]
    required_envs = []
    code_input_sig = [('smiles1', 'str', 'N/A', 'SMILES string of the first molecule'), ('smiles2', 'str', 'N/A', 'SMILES string of the second molecule.')]
    text_input_sig = [('smiles_pair', 'str', 'N/A', 'SMILES strings of the two molecules, separated by a semicolon.')]
    output_sig = [('similarity', 'str', 'Tanimoto similarity score and similarity description')]
    examples = [
        {'code_input': {'smiles1': 'CCO', 'smiles2': 'CCN'}, 'text_input': {'smiles_pair': 'CCO;CCN'}, 'output': {'similarity': 'The Tanimoto similarity between CCO and CCN is 0.3333, indicating that the two molecules are not similar.'}},
        {'code_input': {'smiles1': 'CCO', 'smiles2': 'C(O)C'}, 'text_input': {'smiles_pair': 'CCO;C(O)C'}, 'output': {'similarity': 'Input Molecules Are Identical'}},
    ]
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT")
    ]
    services_and_software = []

    def _run_base(self, smiles1: str, smiles2: str) -> str:
        if not is_smiles(smiles1):
            raise ChemMCPInputError(f"smiles1 `{smiles1}` is not a valid SMILES string.")
        
        if not is_smiles(smiles2):
            raise ChemMCPInputError(f"smiles2 `{smiles2}` is not a valid SMILES string.")

        similarity = tanimoto(smiles1, smiles2)

        if isinstance(similarity, str):
            return similarity

        sim_score = {
            0.9: "very similar",
            0.8: "similar",
            0.7: "somewhat similar",
            0.6: "not very similar",
            0: "not similar",
        }
        if similarity == 1:
            return "The input molecules are identical."
        else:
            val = sim_score[
                max(key for key in sim_score.keys() if key <= round(similarity, 1))
            ]
            message = f"The Tanimoto similarity between {smiles1} and {smiles2} is {round(similarity, 4)}, indicating that the two molecules are {val}."
        return message
    
    def _run_text(self, smiles_pair: str) -> str:
        smiles1, smiles2 = smiles_pair.split(';')
        smiles1 = smiles1.strip()
        smiles2 = smiles2.strip()
        return self._run_base(smiles1, smiles2)


if __name__ == "__main__":
    run_mcp_server()

