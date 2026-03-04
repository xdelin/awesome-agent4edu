import importlib

_tool_module_map = {
    "WebSearch":    "web_search",
    "MoleculeCaptioner": "molecule_captioner",
    "MoleculeGenerator": "molecule_generator",
    "PubchemSearchQA": "pubchem_search_qa",
    "PubchemSearch": "pubchem_search",  # Not registerred as an MCP tool
    "ForwardSynthesis": "forward_synthesis",
    "Retrosynthesis": "retrosynthesis",
    "MoleculeSimilarity": "molecule_similarity",
    "MoleculeWeight": "molecule_weight",
    "FunctionalGroups": "functional_groups",
    "SmilesCanonicalization": "smiles_canonicalization",
    "MoleculeAtomCount": "molecule_atom_count",
    "MoleculePrice": "molecule_price",
    "PatentCheck": "patent_check",
    "SolubilityPredictor": "solubility_predictor",
    "LogDPredictor": "logd_predictor",
    "BbbpPredictor": "bbbp_predictor",
    "ToxicityPredictor": "toxicity_predictor",
    "HivInhibitorPredictor": "hiv_inhibitor_predictor",
    "SideEffectPredictor": "side_effect_predictor",
    "Iupac2Smiles": "iupac2smiles",
    "Smiles2Iupac": "smiles2iupac",
    "Smiles2Formula": "smiles2formula",
    "Name2Smiles": "name2smiles",
    "Selfies2Smiles": "selfies2smiles",
    "Smiles2Selfies": "smiles2selfies",
    "MoleculeSmilesCheck": "molecule_smiles_check",
    "ReactionSmilesCheck": "reaction_smiles_check",
    "Smiles2Cas": "smiles2cas",
    "SafetyCheck": "safety_check",
    "MoleculeModifier": "molecule_modifier",
    "PythonExecutor": "python_executor",
    "MoleculeVisualizer": "molecule_visualizer",
}

__all__ = list(_tool_module_map.keys())

def __getattr__(name: str):
    if name in _tool_module_map:
        module_name = _tool_module_map.get(name)
        if module_name is None:
            raise AttributeError(f"No mapping for tool {name!r} in chemmcp.tools")
        module = importlib.import_module(f"{__name__}.{module_name}")
        try:
            return getattr(module, name)
        except AttributeError:
            raise ImportError(f"Module {module_name!r} has no attribute {name!r}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return list(__all__)
