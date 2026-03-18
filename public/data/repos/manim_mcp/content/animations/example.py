from manim import *

class ExampleScene(Scene):
    def construct(self):
        # Create a circle and a square
        circle = Circle()
        circle.set_fill(BLUE, opacity=0.5)
        circle.set_stroke(BLUE_E, width=4)
        
        square = Square()
        square.set_fill(RED, opacity=0.5)
        square.set_stroke(RED_E, width=4)
        
        # Position the shapes
        VGroup(circle, square).arrange(RIGHT, buff=1)
        
        # Display the shapes
        self.play(Create(circle), Create(square))
        self.wait()
        
        # Animate the shapes
        self.play(
            circle.animate.shift(DOWN),
            square.animate.shift(UP)
        )
        self.wait()
        
        # Transform the square into a star
        star = Star(outer_radius=1, n=5)
        star.set_fill(YELLOW, opacity=0.5)
        star.set_stroke(YELLOW_E, width=4)
        
        self.play(Transform(square, star))
        self.wait()
        
        # Add a title
        title = Text("Manim in Docker").scale(0.75)
        title.to_edge(UP)
        
        self.play(Write(title))
        self.wait(2)


class MathExample(Scene):
    def construct(self):
        # A more advanced example with mathematics
        formula1 = MathTex(r"\frac{d}{dx} \left( \int_{a}^{x} f(t) \, dt \right) = f(x)")
        formula2 = MathTex(r"\frac{1}{2\pi i} \oint_{\gamma} \frac{f(z)}{(z - z_0)^{n+1}}~dz = \frac{f^{(n)}(z_0)}{n!}")
        
        self.play(Write(formula1))
        self.wait(2)
        self.play(ReplacementTransform(formula1, formula2))
        self.wait(2) 