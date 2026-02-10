"""Model validation for Stella system dynamics models."""

import re
from dataclasses import dataclass
from typing import Optional

from .xmile import StellaModel


@dataclass
class ValidationError:
    """Represents a validation error or warning."""
    severity: str  # "error" or "warning"
    category: str  # e.g., "undefined_variable", "mass_balance", "missing_connection"
    message: str
    variable: Optional[str] = None


class ModelValidator:
    """Validates Stella models for common errors."""

    def __init__(self, model: StellaModel):
        self.model = model
        self.errors: list[ValidationError] = []

    def validate(self) -> list[ValidationError]:
        """Run all validation checks and return errors/warnings."""
        self.errors = []

        self._check_undefined_variables()
        self._check_mass_balance()
        self._check_missing_connections()
        self._check_orphan_flows()
        self._check_stock_inflow_outflow_consistency()
        self._check_circular_dependencies()

        return self.errors

    def _get_all_variable_names(self) -> set[str]:
        """Get all variable names in the model."""
        names = set()
        for name in self.model.stocks:
            names.add(name)
        for name in self.model.flows:
            names.add(name)
        for name in self.model.auxs:
            names.add(name)
        return names

    def _extract_variable_references(self, equation: str) -> set[str]:
        """Extract variable names referenced in an equation."""
        if not equation:
            return set()

        # Remove string literals
        equation = re.sub(r'"[^"]*"', '', equation)

        # Remove function names (followed by parentheses)
        # Common Stella functions
        functions = [
            'INIT', 'TIME', 'DT', 'STARTTIME', 'STOPTIME',
            'ABS', 'MIN', 'MAX', 'SUM', 'MEAN', 'SQRT',
            'EXP', 'LN', 'LOG10', 'SIN', 'COS', 'TAN',
            'IF', 'THEN', 'ELSE', 'AND', 'OR', 'NOT',
            'DELAY', 'DELAY1', 'DELAY3', 'SMTH1', 'SMTH3',
            'PULSE', 'STEP', 'RAMP', 'RANDOM', 'NORMAL',
            'GRAPH', 'LOOKUP', 'INTERPOLATE', 'HISTORY',
            'SAFEDIV', 'FORCST', 'TREND', 'NPV', 'IRR',
            'ROUND', 'INT', 'MOD', 'COUNTER', 'PREVIOUS'
        ]

        for func in functions:
            equation = re.sub(rf'\b{func}\s*\(', '(', equation, flags=re.IGNORECASE)

        # Extract potential variable names (alphanumeric with underscores)
        # Match words that could be variable names
        potential_vars = re.findall(r'\b([A-Za-z_][A-Za-z0-9_]*)\b', equation)

        # Filter out numbers and common keywords
        keywords = {'true', 'false', 'pi', 'e', 'inf', 'nan'}
        refs = set()
        for var in potential_vars:
            if var.lower() not in keywords:
                # Check if it's not a pure number
                try:
                    float(var)
                except ValueError:
                    refs.add(var)

        return refs

    def _check_undefined_variables(self):
        """Check for references to undefined variables."""
        all_vars = self._get_all_variable_names()

        # Check flow equations
        for name, flow in self.model.flows.items():
            refs = self._extract_variable_references(flow.equation)
            for ref in refs:
                if ref not in all_vars:
                    self.errors.append(ValidationError(
                        severity="error",
                        category="undefined_variable",
                        message=f"Flow '{flow.name}' references undefined variable '{ref}'",
                        variable=name
                    ))

        # Check aux equations
        for name, aux in self.model.auxs.items():
            refs = self._extract_variable_references(aux.equation)
            for ref in refs:
                if ref not in all_vars:
                    self.errors.append(ValidationError(
                        severity="error",
                        category="undefined_variable",
                        message=f"Auxiliary '{aux.name}' references undefined variable '{ref}'",
                        variable=name
                    ))

    def _check_mass_balance(self):
        """Check for potential mass balance issues."""
        # For each stock, check if it has at least one inflow or outflow
        for name, stock in self.model.stocks.items():
            if not stock.inflows and not stock.outflows:
                self.errors.append(ValidationError(
                    severity="warning",
                    category="mass_balance",
                    message=f"Stock '{stock.name}' has no inflows or outflows (isolated reservoir)",
                    variable=name
                ))

        # Check if any flows reference stocks that don't exist
        for name, flow in self.model.flows.items():
            if flow.from_stock and flow.from_stock not in self.model.stocks:
                self.errors.append(ValidationError(
                    severity="error",
                    category="mass_balance",
                    message=f"Flow '{flow.name}' references non-existent from_stock '{flow.from_stock}'",
                    variable=name
                ))
            if flow.to_stock and flow.to_stock not in self.model.stocks:
                self.errors.append(ValidationError(
                    severity="error",
                    category="mass_balance",
                    message=f"Flow '{flow.name}' references non-existent to_stock '{flow.to_stock}'",
                    variable=name
                ))

    def _check_missing_connections(self):
        """Check for missing connectors based on equation references."""
        all_vars = self._get_all_variable_names()

        # Build set of existing connections
        existing_connections = set()
        for conn in self.model.connectors:
            existing_connections.add((conn.from_var, conn.to_var))

        # Check flow equations for missing connectors
        for name, flow in self.model.flows.items():
            refs = self._extract_variable_references(flow.equation)
            for ref in refs:
                if ref in all_vars and ref != name:
                    if (ref, name) not in existing_connections:
                        self.errors.append(ValidationError(
                            severity="warning",
                            category="missing_connection",
                            message=f"Flow '{flow.name}' uses '{ref}' but no connector exists",
                            variable=name
                        ))

        # Check aux equations for missing connectors
        for name, aux in self.model.auxs.items():
            refs = self._extract_variable_references(aux.equation)
            for ref in refs:
                if ref in all_vars and ref != name:
                    if (ref, name) not in existing_connections:
                        self.errors.append(ValidationError(
                            severity="warning",
                            category="missing_connection",
                            message=f"Auxiliary '{aux.name}' uses '{ref}' but no connector exists",
                            variable=name
                        ))

    def _check_orphan_flows(self):
        """Check for flows that aren't connected to any stock."""
        for name, flow in self.model.flows.items():
            if not flow.from_stock and not flow.to_stock:
                self.errors.append(ValidationError(
                    severity="warning",
                    category="orphan_flow",
                    message=f"Flow '{flow.name}' is not connected to any stock",
                    variable=name
                ))

    def _check_stock_inflow_outflow_consistency(self):
        """Check that stock inflows/outflows match flow definitions."""
        for name, stock in self.model.stocks.items():
            # Check inflows
            for inflow in stock.inflows:
                if inflow not in self.model.flows:
                    self.errors.append(ValidationError(
                        severity="error",
                        category="undefined_variable",
                        message=f"Stock '{stock.name}' references undefined inflow '{inflow}'",
                        variable=name
                    ))
                elif self.model.flows[inflow].to_stock != name:
                    self.errors.append(ValidationError(
                        severity="warning",
                        category="inconsistent_flow",
                        message=f"Stock '{stock.name}' lists '{inflow}' as inflow, but flow doesn't point to this stock",
                        variable=name
                    ))

            # Check outflows
            for outflow in stock.outflows:
                if outflow not in self.model.flows:
                    self.errors.append(ValidationError(
                        severity="error",
                        category="undefined_variable",
                        message=f"Stock '{stock.name}' references undefined outflow '{outflow}'",
                        variable=name
                    ))
                elif self.model.flows[outflow].from_stock != name:
                    self.errors.append(ValidationError(
                        severity="warning",
                        category="inconsistent_flow",
                        message=f"Stock '{stock.name}' lists '{outflow}' as outflow, but flow doesn't originate from this stock",
                        variable=name
                    ))

    def _check_circular_dependencies(self):
        """Check for circular dependencies in auxiliary variables (excluding stocks/flows)."""
        # Build dependency graph for aux variables only
        deps: dict[str, set[str]] = {}
        aux_names = set(self.model.auxs.keys())

        for name, aux in self.model.auxs.items():
            refs = self._extract_variable_references(aux.equation)
            # Only track dependencies on other aux variables
            deps[name] = refs & aux_names

        # Check for cycles using DFS
        def has_cycle(node: str, visited: set[str], rec_stack: set[str]) -> list[str]:
            visited.add(node)
            rec_stack.add(node)

            for neighbor in deps.get(node, set()):
                if neighbor not in visited:
                    cycle = has_cycle(neighbor, visited, rec_stack)
                    if cycle:
                        return [node] + cycle
                elif neighbor in rec_stack:
                    return [node, neighbor]

            rec_stack.remove(node)
            return []

        visited: set[str] = set()
        for name in aux_names:
            if name not in visited:
                cycle = has_cycle(name, visited, set())
                if cycle:
                    cycle_str = " -> ".join(cycle)
                    self.errors.append(ValidationError(
                        severity="error",
                        category="circular_dependency",
                        message=f"Circular dependency detected among auxiliaries: {cycle_str}",
                        variable=cycle[0]
                    ))
                    break  # Only report first cycle found


def validate_model(model: StellaModel) -> list[ValidationError]:
    """Convenience function to validate a model."""
    validator = ModelValidator(model)
    return validator.validate()
