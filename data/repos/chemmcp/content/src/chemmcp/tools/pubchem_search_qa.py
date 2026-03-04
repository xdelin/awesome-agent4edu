from typing import Literal

from ..utils.base_tool import BaseTool
from ..utils.errors import *
from ..tool_utils.llm import llm_completion
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from .pubchem_search import PubchemSearch


QA_SYSTEM_PROMPT = "You are an expert chemist. You will be given the PubChem page about a molecule/compound, and your task is to answer the question based on the information of the page. Your answer should be accurate and concise, and contain all the information necessary to answer the question."



@ChemMCPManager.register_tool
class PubchemSearchQA(BaseTool):
    __version__ = "0.1.1"
    name = "PubchemSearchQA"
    func_name = 'search_pubchem_qa'
    description = "Answer questions about molecules/compounds based on the information from PubChem, one of the most comprehensive database of chemical molecules and their activities. You can get authoritative answers about molecular names, properties, activities, and more."
    implementation_description = "Uses the PubchemSearch tool to get the document of a molecule/compound from [PubChem](https://pubchem.ncbi.nlm.nih.gov/), and prompts an LLM to answer the related question."
    oss_dependencies = []
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/"), ("Custom LLMs", None)]
    categories = ["Molecule"]
    tags = ["PubChem", "Molecular Information", "Molecular Properties", "SMILES", "IUPAC", "Molecular Names", "APIs", "QA", "LLMs"]
    required_envs = [("__llms__", "LiteLLM Envs.")]
    text_input_sig = [("representation_name_and_representation_and_question", "str", 'N/A', "The representation name and representation of the molecule/compound, e.g., \"SMILES: <SMILES>\", \"IUPAC: <IUPAC name>\", or \"Name: <common name>\". Followed by \"Question: <your question about the molecule/compound>\".")]
    code_input_sig = [("representation_name", "str", 'N/A', "The representation name, can be \"smiles\", \"iupac\", or \"name\" (chemical's common name)."), ("representation", "str", 'N/A', "The representation of the molecule/compound, corresponding to the representation_name used."), ("question", "str", 'N/A', "The question about the molecule/compound.")]
    output_sig = [("answer", "str", "The answer to the question based on the PubChem page.")]
    examples = [
        {'code_input': {'representation_name': 'SMILES', 'representation': 'CCO', "question": "What properties do this molecule have?"}, 'text_input': {'representation_name_and_representation_and_question': 'SMILES: CCO   Questions: What properties do this molecule have?'}, 'output': {'answer': 'This molecule has the following properties: [...]'}},
        
        {'code_input': {'representation_name': 'IUPAC', 'representation': 'ethanol', "question": "What properties do this molecule have?"}, 'text_input': {'representation_name_and_representation_and_question': 'IUPAC: ethanol   Questions: What properties do this molecule have?'}, 'output': {'answer': 'This molecule has the following properties: [...]'}},

        {'code_input': {'representation_name': 'Name', 'representation': 'alcohol', "question": "What properties do this molecule have?"}, 'text_input': {'representation_name_and_representation_and_question': 'Name: alcohol   Questions: What properties do this molecule have?'}, 'output': {'answer': 'This molecule has the following properties: [...]'}},
    ]

    def __init__(self, llm_model=None, init=True, interface='text') -> None:
        super().__init__(init, interface)
        self.llm_model = llm_model
        self.pubchem_search = PubchemSearch(init=init, interface='code')

    def _run_text(self, representation_name_and_representation_and_question: str) -> str:
        if 'Question:' not in representation_name_and_representation_and_question:
            raise ChemMCPInputError("The input is not in a correct format. Please input the molecule/compound representation followed by the question about the molecule/compound. An example: \"SMILES: <SMILES of the molecule/compound> Question: <your question about the molecule/compound>\".")
        representation_name_and_representation_and_question, question = representation_name_and_representation_and_question.split('Question:')
        representation_name_and_representation_and_question = representation_name_and_representation_and_question.strip()
        question = question.strip()

        try:
            namespace, identifier = representation_name_and_representation_and_question.split(':')
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
            raise ChemMCPInputError("The input is not in a correct format: %s If searching with SMILES, please input \"SMILES: <SMILES of the molecule/compound>\"; if searching with IUPAC name, please input \"IUPAC: <IUPAC name of the molecule/compound>\"; if searching with common name, please input \"Name: <common name of the molecule/compound>\". After that, append your question about the molecule/compound as \"Question: <your question>\"." % str(e))
        r = self._run_base(namespace, identifier, question)
        return r
    
    def _run_base(self, representation_name: Literal["SMILES", "IUPAC", "Name"], representation: str, question: str):
        doc = self.pubchem_search.run_code(representation_name, representation)
        conversation = [
            {'role': 'system', 'content': QA_SYSTEM_PROMPT},
            {'role': 'user', 'content': doc + '\n\n\n\nQuestion: ' + question},
        ]
        r = llm_completion(messages=conversation, llm_model=self.llm_model)
        r = r.choices[0].message.content.strip()
        return r


if __name__ == "__main__":
    run_mcp_server()
