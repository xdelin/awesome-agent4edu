import logging

from ..utils.errors import ChemMCPToolProcessError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..tool_utils.rxn4chem import RXN4Chem
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ForwardSynthesis(RXN4Chem):
    """Predict reaction."""
    __version__ = "0.1.0"
    name = "ForwardSynthesis"
    func_name = "do_forward_synthesis"
    description = "Given reactants and reagents, predict the product(s) of a chemical reaction."
    implementation_description = "Uses the [IBM RXN for Chem](https://rxn.app.accelerate.science/) API to predict the product(s) of a chemical reaction. Outputs the top 1 prediction."
    categories = ["Reaction"]
    tags = ["Reaction Prediction", "Neural Networks", "APIs", "SMILES"]
    required_envs = [("RXN4CHEM_API_KEY", "The API key for IBM RXN4Chem.")]
    text_input_sig = [("reactants_and_reagents_smiles", "str", "N/A", "The SMILES of the reactants and reagents separated by a dot '.'.")]
    code_input_sig = [("reactants_and_reagents_smiles", "str", "N/A", "The SMILES of the reactants and reagents separated by a dot '.'.")]
    output_sig = [("product_smiles", "str", "The SMILES of the product(s).")]
    examples = [
        {'code_input': {'reactants_and_reagents_smiles': 'CCN.CN1C=CC=C1C=O'}, 'text_input': {'reactants_and_reagents_smiles': 'CCN.CN1C=CC=C1C=O'}, 'output': {'product_smiles': 'CCNCc1cccn1C'}},
    ]
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT"),
        ("rxn4chemistry", "https://github.com/rxn4chemistry/rxn4chemistry", "MIT")
    ]
    services_and_software = [("IBM RXN for Chemistry", "https://rxn.app.accelerate.science/")]

    def _run_text(self, reactants_and_reagents_smiles: str) -> str:
        return self._run_base(reactants_and_reagents_smiles)

    def _run_base(self, reactants_and_reagents_smiles: str) -> str:
        """Run reaction prediction."""
        # Check that input is smiles
        if not is_smiles(reactants_and_reagents_smiles):
            raise ChemMCPInputError("The input contains invalid SMILES. Please double-check.")
        if '.' not in reactants_and_reagents_smiles:
            raise ChemMCPInputError("This tool only support inputs with at least two reactants and reagents separated by a dot '.'. Please double-check.")

        try:
            prediction_id = self.predict_reaction(reactants_and_reagents_smiles)
            results = self.get_results(prediction_id)
            product = results["productMolecule"]["smiles"]
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMCPToolProcessError("Failed to predict the products for the input string. Please make sure the input is a valid SMILES string containing reactants and reagents separated by a dot '.'") from e
        return product

    @RXN4Chem.retry(10, ChemMCPToolProcessError)
    def predict_reaction(self, reactants: str) -> str:
        """Make api request."""
        response = self.rxn4chem.predict_reaction(reactants)
        if "prediction_id" in response.keys():
            return response["prediction_id"]
        else:
            raise ChemMCPToolProcessError("The tool failed to predict the reaction. Maybe the input is invalid. Please make sure the input is valid SMILES of reactants separated by dot '.' and try again.")

    @RXN4Chem.retry(10, ChemMCPToolProcessError)
    def get_results(self, prediction_id: str) -> str:
        """Make api request."""
        results = self.rxn4chem.get_predict_reaction_results(prediction_id)
        if "payload" in results["response"].keys():
            return results["response"]["payload"]["attempts"][0]
        else:
            raise ChemMCPToolProcessError("Error in obtaining the results. Please make sure the input is valid SMILES of reactants separated by dot '.' and try again.")
        

if __name__ == "__main__":
    run_mcp_server()
