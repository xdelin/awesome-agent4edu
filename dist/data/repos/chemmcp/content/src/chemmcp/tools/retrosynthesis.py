import logging
from time import sleep

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPToolProcessError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..tool_utils.rxn4chem import RXN4Chem
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)



@ChemMCPManager.register_tool
class Retrosynthesis(RXN4Chem):
    """Predict single-step retrosynthesis."""
    __version__ = "0.1.0"
    name = "Retrosynthesis"
    func_name = "do_retrosynthesis"
    description = "Conduct single-step retrosynthesis. Given the product(s), predict multiple sets of potential reactants, along with their confidence."
    implementation_description = "Uses the [IBM RXN for Chem](https://rxn.app.accelerate.science/) API to do single-step retrosynthesis. Outputs all possible reactants and their confidence, sorted by confidence in descending order."
    categories = ["Reaction"]
    tags =["Reaction Prediction", "Neural Networks", "APIs", "SMILES"]
    required_envs = [("RXN4CHEM_API_KEY", "The API key for IBM RXN4Chem.")]
    code_input_sig = [("product_smiles", "str", "N/A", "The SMILES of the product.")]
    text_input_sig = [("product_smiles", "str", "N/A", "The SMILES of the product.")]
    output_sig = [("reactants_and_confidence", "str", "The SMILES of the reactants and the confidence.")]
    examples = [
        {'code_input': {'product_smiles': 'CCO'}, 'text_input': {'product_smiles': 'CCO'}, 'output': {'reactants_and_confidence': 'There are 13 possible sets of reactants for the given product:\n1.\tReactants: C1CCOC1.CCNC(=O)c1cccn1C.[Li][AlH4]\tConfidence: 1.0\n2.\tReactants: CCN.CCO.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n3.\tReactants: CCN.CO.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n4.\tReactants: CCN.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n5.\tReactants: CCN.CCO.Cn1cccc1C=O.O.[BH4-].[Na+]\tConfidence: 1.0\n6.\tReactants: CCN.CO.Cn1cccc1C=O.O.[BH4-].[Na+]\tConfidence: 1.0\n7.\tReactants: C1CCOC1.CCN.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n8.\tReactants: CCN.Cl.Cn1cccc1C=O\tConfidence: 0.938\n9.\tReactants: CCN.Cn1cccc1C=O\tConfidence: 0.917\n10.\tReactants: CCN.Cl.Cn1cccc1C=O\tConfidence: 0.841\n11.\tReactants: C1CCOC1.CCN.Cn1cccc1C=O\tConfidence: 0.797\n12.\tReactants: C1CCOC1.CCN.CO.Cn1cccc1C=O\tConfidence: 0.647\n13.\tReactants: C1CCOC1.CC(=O)NCc1cccn1C.[Li][AlH4]\tConfidence: 1.0\n'}},  
    ]
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT"),
        ("rxn4chemistry", "https://github.com/rxn4chemistry/rxn4chemistry", "MIT")
    ]
    services_and_software = [("IBM RXN for Chemistry", "https://rxn.app.accelerate.science/")]

    def _run_base(self, product_smiles: str) -> str:
        """Run retrosynthesis prediction."""
        # Check that input is smiles
        if not is_smiles(product_smiles):
            raise ChemMCPInputError("The input contains invalid SMILES. Please double-check.")

        try:
            prediction_id = self.predict_retrosynthesis(product_smiles)
            paths = self.get_paths(prediction_id)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMCPToolProcessError("Failed to predict the reactants for the input string. Please make sure the input is a valid SMILES string containing one chemical.") from e

        result = "There %s %d possible sets of reactants for the given product:\n" % (
            "are" if len(paths) > 1 else "is",
            len(paths),
        )
        result_list = []
        for idx, path in enumerate(paths, start=1):
            children_smiles, confidence = self._get_children_smiles_and_confidence(path)
            result_list.append((children_smiles, confidence))
        result_list.sort(key=lambda x: x[1], reverse=True)
        for idx, (children_smiles, confidence) in enumerate(result_list, start=1):
            result += f"{idx}.\tReactants: {children_smiles}\tConfidence: {confidence}\n"
        return result

    @RXN4Chem.retry(10, KeyError)
    def predict_retrosynthesis(self, target: str) -> str:
        """Make api request."""
        response = self.rxn4chem.predict_automatic_retrosynthesis(
            product=target,
            max_steps=1,
        )
        if "prediction_id" in response.keys():
            return response["prediction_id"]
        raise KeyError

    @RXN4Chem.retry(20, ChemMCPToolProcessError)
    def get_paths(self, prediction_id: str) -> str:
        """Make api request."""
        results = self.rxn4chem.get_predict_automatic_retrosynthesis_results(
            prediction_id
        )
        if "retrosynthetic_paths" not in results.keys():
            raise ChemMCPToolProcessError("Error in obtaining the results. Please make sure the input is valid SMILES and try again.")
        paths = results["retrosynthetic_paths"]
        if paths is not None:
            if len(paths) > 0:
                return paths
        if results["status"] == "PROCESSING":
            sleep(self.sleep_time * 2)
        raise ChemMCPToolProcessError("Error in obtaining the results. Please make sure the input is valid SMILES and try again.")
    
    def _get_children_smiles_and_confidence(self, path):
        children = path['children']
        children_smiles = []
        for child in children:
            smiles = child['smiles']
            children_smiles.append(smiles)
        return '.'.join(children_smiles), path['confidence']


if __name__ == "__main__":
    run_mcp_server()
