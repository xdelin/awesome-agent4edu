import logging

from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.property_prediction import PropertyPredictor


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SideEffectPredictor(PropertyPredictor):
    subtask_list = ['Blood and lymphatic system disorders', 'Cardiac disorders', 'Congenital, familial and genetic disorders', 'Ear and labyrinth disorders', 'Endocrine disorders', 'Eye disorders', 'Gastrointestinal disorders', 'Hepatobiliary disorders', 'Immune system disorders', 'Metabolism and nutrition disorders', 'Musculoskeletal and connective tissue disorders', 'Neoplasms benign, malignant and unspecified (incl cysts and polyps)', 'Nervous system disorders', 'Pregnancy, puerperium and perinatal conditions', 'Psychiatric disorders', 'Renal and urinary disorders', 'Reproductive system and breast disorders', 'Respiratory, thoracic and mediastinal disorders', 'Skin and subcutaneous tissue disorders', 'Vascular disorders']

    __version__ = "0.1.0"
    name = "SideEffectPredictor"
    func_name = 'predict_side_effect'
    description = "Predict whether a molecule can cause 20 different side effects, along with the probabilities of each side effect. The side effects are: " + '; '.join(['(%d) %s' % (idx, subtask) for idx, subtask in enumerate(subtask_list, start=1)]) + "."
    implementation_description = "Uses the [Uni-Mol](https://github.com/deepmodeling/Uni-Mol) model fine-tuned on [SmolInstruct](https://huggingface.co/datasets/osunlp/SMolInstruct) PP-SIDER data to predict the probabilities of the compound to cause different side effects, and uses a text template to construct textual output."
    categories = ["Molecule"]
    tags = ["Molecular Information", "Molecular Properties", "SMILES", "Neural Networks"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('side_effect', 'str', 'The probabilities of the compound to cause different side effects.')]
    examples = [
        {'code_input': {'smiles': 'CC1=CC(C)=C(NC(=O)CN(CC(=O)O)CC(=O)O)C(C)=C1Br'}, 'text_input': {'smiles': 'CC1=CC(C)=C(NC(=O)CN(CC(=O)O)CC(=O)O)C(C)=C1Br'}, 'output': {'side_effect': "The probabilities of the compound to cause different side effects are as follows:Blood and lymphatic system disorders: 11.29%, which means it's unlikely to cause the side effect.\nCardiac disorders: 10.92%, which means it's unlikely to cause the side effect.\nCongenital, familial and genetic disorders: 11.98%, which means it's unlikely to cause the side effect.\nEar and labyrinth disorders: 8.48%, which means it's unlikely to cause the side effect.\nEndocrine disorders: 4.16%, which means it's unlikely to cause the side effect.\nEye disorders: 15.19%, which means it's unlikely to cause the side effect.\nGastrointestinal disorders: 57.00%, which means it's likely to cause the side effect.\nHepatobiliary disorders: 9.62%, which means it's unlikely to cause the side effect.\nImmune system disorders: 10.14%, which means it's unlikely to cause the side effect.\nMetabolism and nutrition disorders: 15.41%, which means it's unlikely to cause the side effect.\nMusculoskeletal and connective tissue disorders: 10.77%, which means it's unlikely to cause the side effect.\nNeoplasms benign, malignant and unspecified (incl cysts and polyps): 4.92%, which means it's unlikely to cause the side effect.\nNervous system disorders: 34.37%, which means it's unlikely to cause the side effect.\nPregnancy, puerperium and perinatal conditions: 3.32%, which means it's unlikely to cause the side effect.\nPsychiatric disorders: 8.06%, which means it's unlikely to cause the side effect.\nRenal and urinary disorders: 10.64%, which means it's unlikely to cause the side effect.\nReproductive system and breast disorders: 4.59%, which means it's unlikely to cause the side effect.\nRespiratory, thoracic and mediastinal disorders: 16.48%, which means it's unlikely to cause the side effect.\nSkin and subcutaneous tissue disorders: 53.97%, which means it's likely to cause the side effect.\nVascular disorders: 18.45%, which means it's unlikely to cause the side effect.\nNote that the results are predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed."}},
    ]
    oss_dependencies = [
        ("Uni-Mol", "https://github.com/deepmodeling/Uni-Mol", "MIT"),
        ("Uni-Core", "https://github.com/dptech-corp/Uni-Core", "MIT")
    ]
    services_and_software = []

    def __init__(
        self, 
        init=True, 
        interface='text'
    ):
        super().__init__('sider', init, interface=interface)

    def _run_base(self, smiles: str) -> str:
        r = super()._run_base(smiles)
        text = []
        for idx, prob in enumerate(r):
            prob = prob * 100
            text.append(f'{self.subtask_list[idx]}: {prob:.2f}%, which means it\'s {"likely" if prob >= 50 else "unlikely"} to cause the side effect.')
        description = 'The probabilities of the compound to cause different side effects are as follows:\n'
        text = description + '\n'.join(text)
        text += '\nNote that the results are predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'
        return text


if __name__ == "__main__":
    run_mcp_server()
