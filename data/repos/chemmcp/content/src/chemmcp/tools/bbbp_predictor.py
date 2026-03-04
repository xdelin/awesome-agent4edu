import logging

from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.property_prediction import PropertyPredictor


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class BbbpPredictor(PropertyPredictor):
    __version__ = "0.1.0"
    name = "BbbpPredictor"
    func_name = 'predict_bbbp'
    description = 'Predict the blood-brain barrier penetration of a molecule given its SMILES representation.'
    implementation_description = "Use the [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) model fine-tuned on [SmolInstruct](https://huggingface.co/datasets/osunlp/SMolInstruct) PP-BBBP data to predict the BBBP probability, and use a text template to construct textual output."
    categories = ["Molecule"]
    tags = ["Molecular Information", "Molecular Properties", "SMILES", "Neural Networks"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('bbbp', 'str', 'The probability of the compound to penetrate the blood-brain barrier.')]
    examples = [
        {'code_input': {'smiles': 'CCNC(=O)/C=C/C1=CC=CC(Br)=C1'}, 'text_input': {'smiles': 'CCNC(=O)/C=C/C1=CC=CC(Br)=C1'}, 'output': {'bbbp': 'The probability of the compound to penetrate the blood-brain barrier is 99.90%, which means it\'s likely to happen.\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'}},
    ]
    oss_dependencies = [
        ("Uni-Mol", "https://github.com/deepmodeling/Uni-Mol", "MIT"),
        ("Uni-Core", "https://github.com/dptech-corp/Uni-Core", "MIT")
    ]
    services_and_software = []

    def __init__(
        self, 
        init=True, 
        interface='code'
    ):
        super().__init__('bbbp', init, interface=interface)

    def _run_base(self, smiles: str) -> str:
        r = super()._run_base(smiles)
        return 'The probability of the compound to penetrate the blood-brain barrier is {:.2f}%, which means it\'s {} to happen.'.format(r * 100, 'likely' if r >= 0.5 else 'unlikely') + '\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'



if __name__ == "__main__":
    run_mcp_server()
