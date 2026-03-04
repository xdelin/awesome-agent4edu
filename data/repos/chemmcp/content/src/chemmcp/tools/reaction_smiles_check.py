import logging

from ..utils.base_tool import BaseTool
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.smiles import is_smiles


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ReactionSmilesCheck(BaseTool):
    __version__ = "0.1.0"
    name = "ReactionSmilesCheck"
    func_name = 'check_reaction_smiles'
    description = "Check the syntactical validity of a reaction SMILES string ([reactant SMILES]>[reagent SMILES]>[product SMILES])."
    implementation_description = "Validates a reaction SMILES string by checking its format (reactants > reagents > products) and verifying that each component is a valid SMILES string. Provides detailed feedback about which parts of the reaction are invalid if any."
    oss_dependencies = []
    categories = ["Molecule", "Reaction"]
    tags = ["SMILES", "SMARTS", "RDKit", "Molecular Information", "Reaction Information"]
    required_envs = []
    text_input_sig = [("smiles", "str", "N/A", "The SMILES string of a chemical reaction to check.")]
    code_input_sig = [("smiles", "str", "N/A", "The SMILES string of a chemical reaction to check.")]
    output_sig = [("result", "str", "Description of the validity of the SMILES string.")]
    examples = [
        {
            'text_input': {'smiles': 'B.C1=CCCCC1.C1=CCCCC1>>B(C1CCCCC1)C1CCCCC1'}, 
            'code_input': {'smiles': 'B.C1=CCCCC1.C1=CCCCC1>>B(C1CCCCC1)C1CCCCC1'}, 
            'output': {'result': 'The reaction SMILES string is valid.'}
        },
    ]
    services_and_software = []

    def _run_base(self, smiles: str) -> str:
        parts = smiles.split('>')
        if len(parts) != 3:
            return f'The SMILES string contains ">", which indicates a reaction, but it is not a valid reaction because it contains {len(parts)} parts instead of 3 (reactants > reagents > products).'
        
        reactants, reagents, products = parts
        reactants = reactants.strip()
        reagents = reagents.strip()
        products = products.strip()

        reactants_valid = is_smiles(reactants)
        if reagents == '':
            reagents_valid = True
        else:
            reagents_valid = is_smiles(reagents)
        products_valid = is_smiles(products)
        
        if reactants_valid and reagents_valid and products_valid:
            return "The reaction SMILES string is valid."
        else:
            invalid_parts = []
            if not reactants_valid:
                invalid_parts.append('reactant(s)')
            if not reagents_valid:
                invalid_parts.append('reagent(s)')
            if not products_valid:
                invalid_parts.append('product(s)')
            return "The reaction SMILES string is invalid. Specifically, the %s is/are not valid SMILES string(s)." % ', '.join(invalid_parts)


if __name__ == "__main__":
    run_mcp_server()
