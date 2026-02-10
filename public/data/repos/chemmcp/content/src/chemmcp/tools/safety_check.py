import re
import logging
from typing import Optional, Literal

import pandas as pd
import pubchempy as pcp

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPSearchFailError, ChemMCPInputError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.smiles import tanimoto
from ..tool_utils.llm import llm_completion
from .pubchem_search import PubchemSearch
from .iupac2smiles import Iupac2Smiles
from .name2smiles import Name2Smiles
from ..tool_utils.canonicalization import canonicalize_molecule_smiles


logger = logging.getLogger(__name__)


SAFETY_SUMMARY_PROMPT = (
    "You will be given a document that contains information about a chemical compound. Your task is to read the document and write a summary of important health, laboratory, and environemntal safety information."
    'Focus on summarizing the following points, and follow the format "Name: description".'
    "Operator safety: Does this substance represent any danger to the person handling it? What are the risks? What precautions should be taken when handling this substance?"
    "GHS information: What are the GHS signal (hazard level: dangerous, warning, etc.) and GHS classification? What do these GHS classifications mean when dealing with this substance?"
    "Environmental risks: What are the environmental impacts of handling this substance."
    "Societal impact: What are the societal concerns of this substance? For instance, is it a known chemical weapon, is it illegal, or is it a controlled substance for any reason?"
    "For each point, use maximum two sentences. Use only the information provided in the document below."
    "If there is not enough information in a category, you may fill in with your knowledge that you are very confident about, and explicitly state so."
    "\n\n\nHere is the document:\n\n{data}"
)


@ChemMCPManager.register_tool
class SafetyCheck(BaseTool):
    __version__ = "0.1.0"
    name = "SafetyCheck"
    func_name = 'check_safety'
    description = "Check online and generate a summary covering its safety information, including its GHS classification, health and laboratory safety, and societal concerns."
    implementation_description = "Uses PubChem's safety data to provide comprehensive safety information about a molecule, including hazard statements, precautionary statements, and other safety-related data. Also compares the molecule with a database of controlled chemicals to check if it is a controlled chemical or similar to a controlled chemical."
    categories = ["Molecule"]
    tags = ["Safety Checking", "LLMs", "Neural Networks", "APIs", "PubChem"]
    required_envs = [("__llms__", "LiteLLM Envs.")]
    text_input_sig = [("representation_name_and_representation", "str", "N/A", "The representation name and representation of the molecule/compound, e.g., \"SMILES: <SMILES>\", \"IUPAC: <IUPAC name>\", or \"Name: <common name>\".")]
    code_input_sig = [("representation_name", "str", "N/A", "The representation name, can be \"SMILES\", \"IUPAC\", or \"Name\" (chemical's common name)."), ("representation", "str", "N/A", "The representation of the molecule/compound, corresponding to the representation_name used.")]
    output_sig = [("safety_summary", "str", "The safety summary of the molecule/compound.")]
    examples = [
        {'text_input': {'representation_name_and_representation': 'SMILES: CCO'}, 'code_input': {'representation_name': 'SMILES', 'representation': 'CCO'}, 'output': {'safety_summary': '**Name: Operator Safety**: Ethanol is highly flammable and can cause eye and skin irritation upon contact. It is recommended to use personal protective equipment, ensure good ventilation, and avoid any ignition sources when handling ethanol. \n\n**Name: GHS Information**: Ethanol has the GHS signal word "Danger" and is classified as Flammable Liquid Category 2 and Eye Irritation Category 2. These classifications mean that ethanol requires precautions to prevent fires and protect the eyes from irritation, such as using flame-resistant containers and eyewear protection.\n\n**Name: Environmental Risks**: Ethanol is moderately harmful to aquatic life with long-term exposure, suggesting caution to prevent its release into waterways. \n\n**Name: Societal Impact**: Ethanol is not a controlled substance but is widely used in alcoholic beverages, regulated for consumption due to its psychoactive effects. It is not considered a chemical weapon or illegal in terms of chemical safety.'}},
    ]
    oss_dependencies = [
        ("ChemCrow", "https://github.com/ur-whitelab/chemcrow-public", "MIT")
    ]
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/"), ("Custom LLMs", None)]

    def __init__(self, llm_model: Optional[str] = None, init=True, interface='code'):
        super().__init__(init, interface=interface)

        self.llm_model = llm_model

        self.pubchem_search = PubchemSearch()
        self.iupac2smiles = Iupac2Smiles()
        self.name2smiles = Name2Smiles()

        self.cw_df = pd.read_csv('https://raw.githubusercontent.com/OSU-NLP-Group/ChemMCP/refs/heads/main/data/chem_wep_smi_canonicalized.csv')

    def _check_controlled_chemical(self, smiles: str) -> bool:
        smiles = canonicalize_molecule_smiles(
            smiles,
            return_none_for_error=False,
            isomeric=False,
            kekulization=False,
            keep_atom_map=False,
        )
        # Check if the molecule is a controlled chemical
        query_esc = re.escape(smiles)
        found = (
            self.cw_df["canonical_smiles"]
            .astype(str)
            .str.contains(f"^{query_esc}$", regex=True)
            .any()
        )

        return found

    def _check_similar_chemical(self, smiles: str) -> bool:
        max_sim = self.cw_df["canonical_smiles"].apply(lambda x: tanimoto(smiles, x)).max()
        return max_sim > 0.35

    def _run_base(self, namespace: Literal["SMILES", "IUPAC", "Name"], identifier: str) -> str:
        # Get PubChem doc
        try:    
            cid: Optional[int] = self.pubchem_search._search_cid(namespace, identifier)
        except ChemMCPSearchFailError as e:
            logger.debug(f"Looking up PubChem failed: {namespace}: {identifier}")
            cid = None
            pubchem_doc = None
        else:
            pubchem_doc = self.pubchem_search.get_cid_doc_text(cid)

        # Get SMILES
        lower_namespace = namespace.lower()
        if lower_namespace == 'smiles':
            smiles = identifier
        else:
            if cid is not None:
                compound = pcp.Compound.from_cid(cid)
                smiles = compound.isomeric_smiles
            else:
                if lower_namespace == 'iupac':
                    smiles = self.iupac2smiles.run_code(identifier)
                elif lower_namespace == 'name':
                    smiles = self.name2smiles.run_code(identifier)
                else:
                    raise ChemMCPInputError(f"Invalid namespace: `{namespace}`. Please use `SMILES`, `IUPAC`, or `Name`.")
        
        # Check if the molecule is a controlled chemical
        is_controlled = self._check_controlled_chemical(smiles)

        # Check if the molecule is similar to a controlled chemical
        is_similar = self._check_similar_chemical(smiles)

        # Generate a description of the molecule's identify or similarity to a controlled chemical
        if is_controlled:
            description = f"Based on the comparison with an existing database, the compound is a controlled chemical."
        elif is_similar:
            description = f"Based on the comparison with an existing database, the compound is similar to a controlled chemical, with a similarity score > 0.35."
        else:
            description = f"Based on the comparison with an existing database, the compound is not a controlled chemical or similar to a controlled chemical."
        
        if pubchem_doc is None:
            return "Cannot find related information for the given compound from PubChem. We only know that: " + description

        overall_doc = ("" if pubchem_doc is None else pubchem_doc) + '\n\n\n\n\n' + description
        overall_doc = overall_doc.strip()

        conversation = [
            {'role': 'user', 'content': SAFETY_SUMMARY_PROMPT.format(data=overall_doc)},
        ]
        summary = llm_completion(messages=conversation, llm_model=self.llm_model)
        summary = summary.choices[0].message.content.strip()
        
        return summary
    
    def _run_text(self, representation_name_and_representation: str) -> str:
        try:
            namespace, identifier = representation_name_and_representation.split(':')
            namespace = namespace.strip()
            identifier = identifier.strip()
            if namespace.lower() in ('smiles',):
                namespace = 'smiles'
            elif namespace.lower() in ('iupac', 'iupac name'):
                namespace = 'iupac'
            elif namespace.lower() in ('name', 'common name'):
                namespace = 'name'
            elif namespace == '':
                raise ChemMCPInputError('Empty representation name.')
            else:
                raise ChemMCPInputError('The representation name \"%s\" is not supported. Please use \"SMILES\", \"IUPAC\", or \"Name\".' % namespace)
        except (ChemMCPInputError, ValueError) as e:
            raise ChemMCPInputError("The input is not in a correct format: %s If searching with SMILES, please input \"SMILES: <SMILES of the molecule/compound>\"; if searching with IUPAC name, please input \"IUPAC: <IUPAC name of the molecule/compound>\"; if searching with common name, please input \"Name: <common name of the molecule/compound>\"." % str(e))
        r = self._run_base(namespace, identifier)
        return r


if __name__ == "__main__":
    run_mcp_server()
