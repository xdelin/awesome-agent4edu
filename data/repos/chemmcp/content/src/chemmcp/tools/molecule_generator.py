from transformers import T5Tokenizer, T5ForConditionalGeneration

from ..utils.base_tool import BaseTool
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class MoleculeGenerator(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeGenerator"
    func_name = "generate_molecule_from_description"
    description = "Generate a molecule represented in SMILES with MolT5 that matches the given textual description."
    implementation_description = "Uses [the MolT5-large model](laituan245/molt5-large-caption2smiles), a transformer-based neural network trained on molecule-text pairs, to generate SMILES representations from natural language descriptions."
    categories = ["Molecule"]
    tags = ["Molecular Information", "Text", "Neural Networks", "SMILES"]
    required_envs = []
    text_input_sig = [('description', 'str', 'N/A', 'Textual description of the molecule.')]
    code_input_sig = [('description', 'str', 'N/A', 'Textual description of the molecule.')]
    output_sig = [('smiles', 'str', 'SMILES representation of the molecule.')]
    examples = [
        {'code_input': {'description': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.'}, 'text_input': {'description': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.'}, 'output': {'smiles': 'CCO\n\nNote: This is a generated SMILES and may not be accurate. Please double check the result.'}},
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
        tokenizer = T5Tokenizer.from_pretrained("laituan245/molt5-large-caption2smiles", model_max_length=512)
        model = T5ForConditionalGeneration.from_pretrained('laituan245/molt5-large-caption2smiles')
        return tokenizer, model
    
    def _run_molt5(self, text):
        if self.tokenizer is None or self.model is None:
            self._init_modules()
        input_ids = self.tokenizer(text, return_tensors="pt").input_ids
        outputs = self.model.generate(input_ids, num_beams=5, max_length=512)
        smiles = self.tokenizer.decode(outputs[0], skip_special_tokens=True)
        return smiles
    
    def _run_base(self, description: str) -> str:
        return self._run_molt5(description) + "\n\nNote: This is a generated SMILES and may not be accurate. Please double check the result."


if __name__ == "__main__":
    run_mcp_server()
