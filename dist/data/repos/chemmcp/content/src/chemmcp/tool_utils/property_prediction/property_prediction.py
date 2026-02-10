import os

from ...utils.base_tool import BaseTool
from ...utils.errors import ChemMCPInputError, ChemMCPToolInitError
from ...tool_utils.smiles import is_smiles
from . import utils as pp_utils


file_path = os.path.abspath(__file__)
workdir = os.getcwd()
rel_path = os.path.relpath(file_path, workdir)
dir_path = os.path.dirname(rel_path)
if dir_path == '':
    dir_path = '.'


MODEL_ARGS = {
    'esol': {
        'model_loc': os.path.join(dir_path, 'checkpoints/esol/checkpoint_best.pt'),
        'task_name': 'esol',
        'task_num': 1,
        'loss_func': 'finetune_mse',
    },
    'lipo': {
        'model_loc': os.path.join(dir_path, 'checkpoints/lipo/checkpoint_best.pt'),
        'task_name': 'lipo',
        'task_num': 1,
        'loss_func': 'finetune_mse',
    },
    'bbbp': {
        'model_loc': os.path.join(dir_path, 'checkpoints/bbbp/checkpoint_best.pt'),
        'task_name': 'bbbp',
    },
    'clintox': {
        'model_loc': os.path.join(dir_path, 'checkpoints/clintox/checkpoint_best.pt'),
        'task_name': 'clintox',
        'task_num': 1,
        'loss_func': 'multi_task_BCE',
    },
    'hiv': {
        'model_loc': os.path.join(dir_path, 'checkpoints/hiv/checkpoint_best.pt'),
        'task_name': 'hiv',
    },
    'sider': {
        'model_loc': os.path.join(dir_path, 'checkpoints/sider/checkpoint_best.pt'),
        'task_name': 'sider',
        'task_num': 20,
        'loss_func': 'multi_task_BCE',
    },
}


class PropertyPredictor(BaseTool):
    __abstract__ = True


    def __init__(
        self, 
        task_name,
        init=True, 
        interface='text'
    ):
        self.task_name = task_name
        self.model = None
        self.task = None
        self.loss = None
        self.args = None
        super().__init__(init, interface=interface)

    def _init_modules(self):
        task_name = self.task_name
        cmd = pp_utils.construct_cmd(**MODEL_ARGS[task_name])
        self.args = pp_utils.parse_args(cmd)
        try:
            self.model, self.task, self.loss = pp_utils.load_model(self.args)
        except FileNotFoundError as e:
            raise ChemMCPToolInitError(f"Model file not found: {e}")

    def _run_base(self, smiles: str) -> str:
        if self.model is None or self.task is None or self.loss is None or self.args is None:
            self._init_modules()

        if not is_smiles(smiles):
            raise ChemMCPInputError(f"Invalid SMILES: {smiles}")
        task_num = 2
        if 'task_num' in MODEL_ARGS[self.task_name]:
            task_num = MODEL_ARGS[self.task_name]['task_num']
        loss_func = 'finetune_cross_entropy'
        if 'loss_func' in MODEL_ARGS[self.task_name]:
            loss_func = MODEL_ARGS[self.task_name]['loss_func']
        r = pp_utils.run_on_smiles(smiles, self.task_name, self.args, self.task, self.model, self.loss, task_num=task_num, loss_func=loss_func)
        return r
