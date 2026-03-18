"""Generate test .stmx files for visual inspection in Stella Professional."""

from pathlib import Path

from stella_mcp.xmile import StellaModel

OUT = Path(__file__).parent


def linear_chain():
    """A → B → C → D: simple linear chain, 3 flows."""
    m = StellaModel("Linear Chain")
    m.add_stock("Atmosphere", "1000")
    m.add_stock("Vegetation", "500")
    m.add_stock("Soil", "2000")
    m.add_stock("Ocean", "38000")

    m.add_flow("GPP", "Atmosphere * 0.01", from_stock="Atmosphere", to_stock="Vegetation")
    m.add_flow("Litter", "Vegetation * 0.02", from_stock="Vegetation", to_stock="Soil")
    m.add_flow("Runoff", "Soil * 0.005", from_stock="Soil", to_stock="Ocean")

    m.add_aux("photosynthesis_rate", "0.01")
    m.add_connector("photosynthesis_rate", "GPP")
    m.add_connector("Atmosphere", "GPP")

    (OUT / "linear_chain.stmx").write_text(m.to_xml())
    print("  linear_chain.stmx")


def feedback_loop():
    """Predator-prey: two stocks with bidirectional flows."""
    m = StellaModel("Predator Prey")
    m.add_stock("Prey", "100")
    m.add_stock("Predators", "20")

    m.add_flow("prey_births", "Prey * prey_birth_rate", to_stock="Prey")
    m.add_flow("predation", "Prey * Predators * predation_rate",
               from_stock="Prey", to_stock="Predators")
    m.add_flow("predator_deaths", "Predators * death_rate", from_stock="Predators")

    m.add_aux("prey_birth_rate", "0.05")
    m.add_aux("predation_rate", "0.001")
    m.add_aux("death_rate", "0.03")

    m.add_connector("prey_birth_rate", "prey_births")
    m.add_connector("Prey", "prey_births")
    m.add_connector("Prey", "predation")
    m.add_connector("Predators", "predation")
    m.add_connector("predation_rate", "predation")
    m.add_connector("death_rate", "predator_deaths")
    m.add_connector("Predators", "predator_deaths")

    (OUT / "feedback_loop.stmx").write_text(m.to_xml())
    print("  feedback_loop.stmx")


def three_stock_cycle():
    """A → B → C → A: circular feedback loop (triangle)."""
    m = StellaModel("Three Stock Cycle")
    m.add_stock("A", "100")
    m.add_stock("B", "50")
    m.add_stock("C", "25")

    m.add_flow("f_ab", "A * 0.1", from_stock="A", to_stock="B")
    m.add_flow("f_bc", "B * 0.1", from_stock="B", to_stock="C")
    m.add_flow("f_ca", "C * 0.1", from_stock="C", to_stock="A")

    m.add_connector("A", "f_ab")
    m.add_connector("B", "f_bc")
    m.add_connector("C", "f_ca")

    (OUT / "three_stock_cycle.stmx").write_text(m.to_xml())
    print("  three_stock_cycle.stmx")


def dense_model():
    """15-element model: 4 stocks, 5 flows, 6 auxs — lots of connectors."""
    m = StellaModel("Dense Model")
    m.add_stock("Population", "1000")
    m.add_stock("Food", "5000")
    m.add_stock("Pollution", "0")
    m.add_stock("Capital", "100")

    m.add_flow("births", "Population * birth_rate", to_stock="Population")
    m.add_flow("deaths", "Population * death_rate", from_stock="Population")
    m.add_flow("consumption", "Population * food_per_capita",
               from_stock="Food", to_stock="Population")
    m.add_flow("emissions", "Population * emission_rate",
               from_stock="Capital", to_stock="Pollution")
    m.add_flow("investment", "Capital * investment_rate", to_stock="Capital")

    m.add_aux("birth_rate", "0.03")
    m.add_aux("death_rate", "0.01")
    m.add_aux("food_per_capita", "0.5")
    m.add_aux("emission_rate", "0.02")
    m.add_aux("investment_rate", "0.05")
    m.add_aux("carrying_capacity", "10000")

    m.add_connector("birth_rate", "births")
    m.add_connector("Population", "births")
    m.add_connector("death_rate", "deaths")
    m.add_connector("Population", "deaths")
    m.add_connector("food_per_capita", "consumption")
    m.add_connector("Population", "consumption")
    m.add_connector("emission_rate", "emissions")
    m.add_connector("Population", "emissions")
    m.add_connector("investment_rate", "investment")
    m.add_connector("Capital", "investment")
    m.add_connector("carrying_capacity", "births")

    (OUT / "dense_model.stmx").write_text(m.to_xml())
    print("  dense_model.stmx")


def user_pinned():
    """Model with some user-specified positions to verify pinning."""
    m = StellaModel("User Pinned Positions")
    m.add_stock("A", "100", x=200, y=300)
    m.add_stock("B", "50", x=600, y=300)
    m.add_stock("C", "25")  # auto-positioned

    m.add_flow("f1", "A * rate", from_stock="A", to_stock="B")
    m.add_flow("f2", "B * 0.1", from_stock="B", to_stock="C")

    m.add_aux("rate", "0.1")
    m.add_connector("rate", "f1")
    m.add_connector("A", "f1")

    (OUT / "user_pinned.stmx").write_text(m.to_xml())
    print("  user_pinned.stmx")


def two_subsystems():
    """Two independent subsystems that should be visually separated."""
    m = StellaModel("Two Subsystems")

    # Subsystem 1: SIR epidemic model
    m.add_stock("Susceptible", "990")
    m.add_stock("Infected", "10")
    m.add_stock("Recovered", "0")
    m.add_flow("infection", "Susceptible * Infected * beta / 1000",
               from_stock="Susceptible", to_stock="Infected")
    m.add_flow("recovery", "Infected * gamma",
               from_stock="Infected", to_stock="Recovered")
    m.add_aux("beta", "0.3")
    m.add_aux("gamma", "0.1")
    m.add_connector("beta", "infection")
    m.add_connector("Susceptible", "infection")
    m.add_connector("Infected", "infection")
    m.add_connector("gamma", "recovery")
    m.add_connector("Infected", "recovery")

    # Subsystem 2: independent economy (no connections to SIR)
    m.add_stock("GDP", "1000")
    m.add_stock("Debt", "500")
    m.add_flow("growth", "GDP * growth_rate", to_stock="GDP")
    m.add_flow("borrowing", "GDP * 0.05", to_stock="Debt")
    m.add_aux("growth_rate", "0.02")
    m.add_connector("growth_rate", "growth")
    m.add_connector("GDP", "growth")
    m.add_connector("GDP", "borrowing")

    (OUT / "two_subsystems.stmx").write_text(m.to_xml())
    print("  two_subsystems.stmx")


if __name__ == "__main__":
    print(f"Generating test .stmx files in {OUT}/")
    linear_chain()
    feedback_loop()
    three_stock_cycle()
    dense_model()
    user_pinned()
    two_subsystems()
    print("Done.")
