from transformers import T5Tokenizer, T5ForConditionalGeneration

from ..utils.base_tool import BaseTool
from ..tool_utils.smiles import is_smiles
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class MoleculeCaptioner(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeCaptioner"
    func_name = "generate_molecule_caption"
    description = "Generate a textual description of the molecule from its SMILES representation with MolT5. This tool uses neural networks to generate descriptions, which may not be accurate or correct. Please first try other tools that provide accurate and authoritative information, and only use this one as the last resort."
    implementation_description = "Uses [the MolT5-large model](laituan245/molt5-large-smiles2caption), a transformer-based neural network trained on molecule-text pairs, to generate natural language descriptions of molecules from their SMILES representations."
    categories = ["Molecule"]
    tags = ["Molecular Information", "Text", "Neural Networks", "SMILES"]
    required_envs = []
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES representation of the molecule.')]
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES representation of the molecule.')]
    output_sig = [('description', 'str', 'Textual description of the molecule.')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'description': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.\n\nNote: This is a generated description and may not be accurate. Please double check the result.'}},
    ]
    oss_dependencies = [
        ("MolT5", "https://github.com/blender-nlp/MolT5", "BSD 3-Clause")
    ]
    services_and_software = []

    def __init__(self, init=True, interface='text') -> None:
        self.tokenizer, self.model = None, None
        super().__init__(init, interface=interface)

    def _init_modules(self):
        self.tokenizer, self.model = self.__load_molt5()
        
    def __load_molt5(self):
        tokenizer = T5Tokenizer.from_pretrained("laituan245/molt5-large-smiles2caption", model_max_length=1024)
        model = T5ForConditionalGeneration.from_pretrained('laituan245/molt5-large-smiles2caption')
        return tokenizer, model
    
    def _run_molt5(self, smiles):
        if self.tokenizer is None or self.model is None:
            self._init_modules()
        input_ids = self.tokenizer(smiles, return_tensors="pt").input_ids
        outputs = self.model.generate(input_ids, num_beams=5, max_length=1024)
        text = self.tokenizer.decode(outputs[0], skip_special_tokens=True)
        return text
    
    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError("The input is not a valid SMILES string.")
        
        return self._run_molt5(smiles) + "\n\nNote: This is a generated description and may not be accurate. Please double check the result."


if __name__ == "__main__":
    run_mcp_server()
