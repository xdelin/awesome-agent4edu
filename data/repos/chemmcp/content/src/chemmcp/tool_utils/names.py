import logging

from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as molD
import pubchempy as pcp
import requests

from ..utils.errors import ChemMCPSearchFailError
from .pubchem import pubchem_iupac2cid, pubchem_name2cid


logger = logging.getLogger(__name__)


def pubchem_iupac2smiles(
    query: str,
    strict: bool = False,
) -> str:
    cid = pubchem_iupac2cid(query, strict=strict)  # May be a tuple in the case of multiple components in the IUPAC name
    if not isinstance(cid, tuple):
        cid = (cid,)
    smiles = []
    for single_cid in cid:
        c = pcp.Compound.from_cid(single_cid)
        r = c.isomeric_smiles
        smiles.append(r.strip())
    r = '.'.join(smiles)

    return r


def pubchem_name2smiles(
    query: str,
) -> str:
    cid = pubchem_name2cid(query)
    c = pcp.Compound.from_cid(cid)
    r = c.isomeric_smiles

    return r


def pubchem_smiles2cas(smi):
    """This function queries the given molecule smiles and returns cas"""

    url_cid = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/{}/cids/JSON"
    url_data = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON"

    try:
        url_cid = url_cid.format('smiles', smi)
        cid = requests.get(url_cid).json()["IdentifierList"]["CID"][0]
        url_data = url_data.format(cid)
        data = requests.get(url_data).json()
    except (requests.exceptions.RequestException, KeyError):
        raise ChemMCPSearchFailError("No Pubchem entry for the input SMILES.")

    try:
        for section in data["Record"]["Section"]:
            if section.get("TOCHeading") == "Names and Identifiers":
                for subsection in section["Section"]:
                    if subsection.get("TOCHeading") == "Other Identifiers":
                        for subsubsection in subsection["Section"]:
                            if subsubsection.get("TOCHeading") == "CAS":
                                return subsubsection["Information"][0]["Value"][
                                    "StringWithMarkup"
                                ][0]["String"]
    except KeyError:
        raise ChemMCPSearchFailError("Cannot find the CAS for the input SMILES.")

    raise ChemMCPSearchFailError("Cannot find the CAS for the input SMILES.")


def pubchem_smiles2iupac(smi):
    """This function queries the given molecule smiles and returns iupac"""

    c = pcp.get_compounds(smi, 'smiles')
    
    if len(c) == 0:
        parts = smi.split('.')
        if len(parts) > 1:
            parts_iupac = []
            parts_cannot_find = []
            for part in parts:
                try:
                    iupac = pubchem_smiles2iupac(part)
                except ChemMCPSearchFailError:
                    parts_cannot_find.append(part)
                else:
                    parts_iupac.append(iupac)
            if len(parts_cannot_find) > 0:
                raise ChemMCPSearchFailError("Cannot find a matched molecule/compound for the following parts of the input SMILES: %s" % ', '.join(parts_cannot_find))
            else:
                r = ';'.join(parts_iupac)
        else:
            raise ChemMCPSearchFailError("Cannot find a matched molecule/compound. Please check the input SMILES.")
    elif len(c) >= 1:
        if len(c) > 1:
            logger.info("There are more than one molecules/compounds that match the input SMILES. Using the first matched one.")
        c = c[0]
    
        r = c.iupac_name
    
    if r is None or r == 'None':
        raise ChemMCPSearchFailError(f"The PubChem entry (CID: {c.cid}) does not have a valid IUPAC name recorded.")

    return r


def addHs(mol):
    mol = Chem.rdmolops.AddHs(mol, explicitOnly=True)
    return mol


def smiles2formula(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = addHs(mol)
    formula = molD.CalcMolFormula(mol)
    return formula
