from rdkit import Chem

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class MoleculeAtomCount(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeAtomCount"
    func_name = 'count_molecule_atoms'
    description = "Count the number of atoms of each type in a molecule."
    implementation_description = "Uses RDKit to parse the SMILES string and count the occurrences of each atom type in the molecule. Returns a detailed description of the atom counts, excluding hydrogen atoms."
    categories = ["Molecule"]
    tags = ["Molecular Information", "RDKit", "SMILES"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('atom_counts', 'str', 'A description of atom numbers in the molecule.')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'atom_counts': 'There are altogether 3 atoms (omitting hydrogen atoms). The types and corresponding numbers are: {\'C\': 2, \'O\': 1}'}},
    ]
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT")
    ]
    services_and_software = []

    def _run_base(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMCPInputError(f"`{smiles}` is not a valid SMILES string.")
        
        num_atoms = mol.GetNumAtoms()
        # Get the atom types
        atom_types = [atom.GetSymbol() for atom in mol.GetAtoms()]
        # Count the occurrences of each atom type
        atom_type_counts = {atom: atom_types.count(atom) for atom in set(atom_types)}
        
        text = "There are altogether %d atoms (omitting hydrogen atoms). The types and corresponding numbers are: %s" % (num_atoms, str(atom_type_counts))

        return text


if __name__ == "__main__":
    run_mcp_server()

