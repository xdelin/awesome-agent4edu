import requests
from copy import deepcopy
from dataclasses import dataclass
from typing import Optional, Literal

import pubchempy as pcp

from ..utils.base_tool import BaseTool
from ..utils.errors import *
from ..tool_utils.smiles import is_smiles
from ..tool_utils.pubchem import pubchem_iupac2cid, pubchem_name2cid


unuseful_section_names = {
    "Structures": None,
    "Chemical Safety": None,
    "Names and Identifiers": {
        "Other Identifiers": None,
        "Synonyms": None,
        "Create Date": None,
        "Modify Date": None,
    },
    "Chemical and Physical Properties": {
        "SpringerMaterials Properties": None,
    },
    "Spectral Information": None,
    "Related Records": None,
    "Chemical Vendors": None,
    "Drug and Medication Information": {
        "WHO Essential Medicines": None,
        "FDA Approved Drugs": None,
        "FDA Orange Book": None,
        "FDA National Drug Code Directory": None,
        "FDA Green Book": None,
        "Drug Labels": None,
        "Clinical Trials": None,
        # "Therapeutic Uses": None,
        # "Drug Warnings": None,
        # "Drug Idiosyncrasies": None,
        # "Reported Fatal Dose": None,
        # "Maximum Drug Dose": None,
    },
    "Pharmacology and Biochemistry": None,
    "Use and Manufacturing": None,
    "Identification": None,
    "Literature": None,
    "Patents": None,
    "Interactions and Pathways": None,
    "Biological Test Results": None,
    "Classification": None,
    "Taxonomy": None,
}


@dataclass
class Information:
    information_item: dict

    @classmethod
    def construct(cls, data):
        information = cls(data)
        return information

    def generate_text(self, display_controls=None):
        data = self.information_item
        value = data['Value']
        text = ""
        if 'StringWithMarkup' in value:
            strings = value['StringWithMarkup']
            for item in strings:
                tmp_text = item['String']
                tmp_unit = (" " + item['Unit']) if 'Unit' in item else ''
                text += tmp_text + tmp_unit + '\n'
        elif 'Number' in value:
            if 'Name' in value:
                name = value['Name']
                text += name + ": "
            strings = value['Number']
            strings = [str(item) for item in strings]
            text += ', '.join(strings)
            tmp_unit = (" " + value['Unit']) if 'Unit' in value else ''
            text += tmp_unit + '\n'
        
        if text.strip() == "":
            return None
        text = text
        return text


@dataclass
class Section:
    level: int
    title: str
    description: str
    display_controls: dict
    information_list: Optional[list] = None
    subsection_list: Optional[list] = None

    @classmethod
    def construct(cls, data, level=1):
        title = data['TOCHeading']
        description = data['Description'] if 'Description' in data else None
        display_controls = data['DisplayControls'] if 'DisplayControls' in data else None

        section = cls(level=level, title=title, description=description, display_controls=display_controls)

        if 'Information' in data:
            information_list = []
            for information_data in data['Information']:
                information = Information.construct(information_data)
                information_list.append(information)
            section.information_list = information_list
        
        if 'Section' in data:
            subsection_list = []
            for subsection_data in data['Section']:
                subsection = Section.construct(subsection_data, level=level + 1)
                subsection_list.append(subsection)
            section.subsection_list = subsection_list

        return section
        
    def generate_text(self, indices=None) -> str:
        # if self.display_controls is not None and "HideThisSection" in self.display_controls and self.display_controls["HideThisSection"] is True:
        #     return None
        
        title_text = '#' * self.level + ' ' + (('.'.join(indices) + ' ') if len(indices) > 0 else '') + self.title + '\n'
        if self.description is not None:
            title_text += 'Section Description: ' + self.description
        title_text += '\n\n'

        if indices is None:
            indices = tuple()

        content_text = ""

        if self.information_list is not None:
            for information in self.information_list:
                tmp_text = information.generate_text(self.display_controls)
                if tmp_text is not None:
                    content_text += tmp_text

        if self.subsection_list is not None:
            idx = 1
            for subsection in self.subsection_list:
                tmp_text = subsection.generate_text(indices + (str(idx),))
                if tmp_text is not None:
                    idx += 1
                    content_text += tmp_text
        
        if content_text.strip() == "":
            return None

        text = title_text + content_text + '\n\n'

        return text


@dataclass
class PubchemStructuredDoc:
    doc_data = []

    @classmethod
    def construct(cls, sections):
        doc = PubchemStructuredDoc()
        section_list = []
        
        for section_data in sections:
            section = Section.construct(section_data)
            section_list.append(section)
        
        doc.doc_data = section_list

        return doc

    def generate_text(self) -> str:
        text = ""
        idx = 1
        for section in self.doc_data:
            tmp_text = section.generate_text(indices=(str(idx),))
            if tmp_text is not None:
                text += tmp_text
                idx += 1
        return text


class PubchemSearch(BaseTool):
    __version__ = "0.1.1"
    name = "PubchemSearch"
    func_name = 'search_pubchem'
    description = "Search for molecule/compound information on PubChem, one of the most comprehensive database of chemical molecules and their activities. You can get authoritative information about molecular names, properties, activities, and more."
    implementation_description = "Uses the PubChem API to search for molecular information using various identifiers (SMILES, IUPAC name, or common name). Filters out unuseful sections and returns markdown formatted information about the molecule including properties, safety data, and other relevant information."
    oss_dependencies = []
    services_and_software = [("PubChem", "https://pubchem.ncbi.nlm.nih.gov/")]
    categories = ["Molecule"]
    tags = ["PubChem", "Molecular Information", "Molecular Properties", "SMILES", "IUPAC", "Molecular Names", "APIs"]
    required_envs = []
    text_input_sig = [("representation_name_and_representation", "str", "N/A", "The representation name and representation of the molecule/compound, e.g., \"SMILES: <SMILES>\", \"IUPAC: <IUPAC name>\", or \"Name: <common name>\".")]
    code_input_sig = [("representation_name", "str", "N/A", "The representation name, can be \"SMILES\", \"IUPAC\", or \"Name\" (chemical's common name)."), ("representation", "str", "N/A", "The representation of the molecule/compound, corresponding to the representation_name used.")]
    output_sig = [("compound_doc", "str", "The document of the molecule/compound in a markdown format.")]
    examples = [
        {'code_input': {'representation_name': 'SMILES', 'representation': 'CCO'}, 'text_input': {'representation_name_and_representation': 'SMILES: CCO'}, 'output': {'compound_doc': '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'}},
        
        {'code_input': {'representation_name': 'IUPAC', 'representation': 'ethanol'}, 'text_input': {'representation_name_and_representation': 'IUPAC: ethanol'}, 'output': {'compound_doc': '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'}},

        {'code_input': {'representation_name': 'Name', 'representation': 'alcohol'}, 'text_input': {'representation_name_and_representation': 'Name: alcohol'}, 'output': {'compound_doc': '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'}},
    ]
    
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/'

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
    
    def _run_base(self, representation_name: Literal["SMILES", "IUPAC", "Name"], representation: str) -> str:
        cid = self._search_cid(representation_name, representation)
        return self.get_cid_doc_text(cid)
    
    def get_data(self, representation_name: Literal["SMILES", "IUPAC", "Name"], representation: str) -> dict:
        cid = self._search_cid(representation_name, representation)
        return self._get_data(cid)
        
    def get_cid_doc_text(self, cid):
        data = self._get_data(cid)

        try:
            sections = self.remove_unuseful_sections(data['Record']['Section'])
        except KeyError:
            print(data)
            print('cid: ', cid)
            raise

        doc = self.construct_doc_text(sections)
        
        return doc

    def _search_cid(self, namespace, identifier):
        namespace = namespace.lower()
        if namespace == 'smiles' and not is_smiles(identifier):
            raise ChemMCPInputError('The input SMILES is invalid. Please double-check. Note that you should input only one molecule/compound at a time.')
        
        if namespace == 'iupac':
            cid = pubchem_iupac2cid(identifier)
        elif namespace == 'smiles':
            try:
                c = pcp.get_compounds(identifier, namespace=namespace)
            except pcp.BadRequestError:
                raise ChemMCPSearchFailError("Error occurred while searching for the molecule/compound on PubChem. Please try other tools or double check your input.")
            if len(c) >= 1:
                c = c[0]
            else:
                raise ChemMCPSearchFailError("Could not find a matched molecule/compound on PubChem. Please double check your input and search for one molecule/compound at a time, or use its another identifier (e.g., IUPAC name or common name) for the search.")
            cid = c.cid
        else:
            cid = pubchem_name2cid(identifier)
        
        if cid is None:
            raise ChemMCPSearchFailError("Could not find a matched molecule/compound on PubChem. Please double check your input and search for one molecule/compound at a time, or use its another identifier for the search.")

        return cid
    
    @staticmethod
    def _get_data(cid):
        url = PubchemSearch.url.format(cid)
        data = requests.get(url).json()
        return data
    
    @staticmethod
    def construct_doc_text(sections):
        doc = PubchemStructuredDoc.construct(sections)
        text = doc.generate_text()
        return text
    
    @staticmethod
    def remove_unuseful_sections(sections):
        sections = deepcopy(sections)

        new_sections = []
        for section in sections:

            section_title = section['TOCHeading']
            if section_title in unuseful_section_names and unuseful_section_names[section_title] is None:
                continue

            if 'Section' in section:
                subsection_list = section['Section']
                new_subsection_list = []
                for subsection in subsection_list:
                    subsection_title = subsection['TOCHeading']
                    if section_title in unuseful_section_names and subsection_title in unuseful_section_names[section_title] and unuseful_section_names[section_title][subsection_title] is None:
                        continue
                    new_subsection_list.append(subsection)
                
                if len(new_subsection_list) == 0:
                    continue

                section['Section'] = new_subsection_list

            new_sections.append(section)
        
        return new_sections
