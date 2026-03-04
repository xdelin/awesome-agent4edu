#!/usr/bin/env python3
"""
LEAP MCP: Basic Explainer Video Creator
==============================================
A minimal Model Context Protocol server for generating simple explainer videos using Manim.
This is a basic implementation - see the full LEAP platform for advanced features.
"""

import re
import shutil
import subprocess
import tempfile
import warnings
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Tuple

from fastmcp import FastMCP
from pydantic import BaseModel, Field, field_validator

# Suppress non-critical warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)


# CONSTANTS AND CONFIGURATION


class Config:
    """Application configuration constants."""
    
    # Server metadata
    SERVER_NAME = "leap-mcp"
    VERSION = "1.0.0"
    
    # Directory structure
    SCRIPTS_DIR = "generated_scripts"
    MEDIA_DIR = "media/videos"
    
    # Manim settings
    MANIM_QUALITY = "medium"  # low, medium, high, fourk
    MANIM_TIMEOUT = 300  # seconds
    DEFAULT_VOICE = "nova"
    
    # Video structure timing (in seconds)
    SCENE_TIMINGS = {
        "hook": (15, 20),
        "phenomenon": (30, 40),
        "mechanism": (40, 50),
        "synthesis": (20, 30)
    }
    
    # Collision avoidance settings
    COLLISION_BUFFER = 0.1
    SPIRAL_SEARCH_MAX_RADIUS = 3.0
    SPIRAL_SEARCH_STEP = 0.3

class VoiceOption(str, Enum):
    """Available OpenAI TTS voices."""
    ALLOY = "alloy"
    ECHO = "echo"
    FABLE = "fable"
    NOVA = "nova"
    ONYX = "onyx"
    SHIMMER = "shimmer"

# ============================================================================
# DATA MODELS
# ============================================================================

class VideoRequest(BaseModel):
    """Model for video creation requests."""
    
    topic: str = Field(..., description="The educational topic to explain")
    voice: VoiceOption = Field(default=VoiceOption.NOVA, description="OpenAI TTS voice")
    
    @field_validator('topic')
    def validate_topic(cls, v):
        """Ensure topic is not empty and has reasonable length."""
        if not v or not v.strip():
            raise ValueError("Topic cannot be empty")
        if len(v) > 200:
            raise ValueError("Topic too long (max 200 characters)")
        return v.strip()

class ExecutionRequest(BaseModel):
    """Model for video execution requests."""
    
    scene_code: str = Field(..., description="Complete Manim scene code")
    topic: str = Field(default="educational", description="Topic for naming")
    
    @field_validator('scene_code')
    def validate_code(cls, v):
        """Basic validation of scene code."""
        if not v or len(v.strip()) < 50:
            raise ValueError("Scene code appears to be too short or empty")
        return v

# COLLISION AVOIDANCE SYSTEM


class CollisionAvoidanceCode:
    """Manages the collision avoidance code injection."""
    
    @staticmethod
    def get_code() -> str:
        """Return optimized collision avoidance code for injection."""
        return '''# Collision Avoidance System (Auto-injected by LEAP MCP)
import numpy as np
from typing import List, Optional, Tuple

class CollisionAvoidanceManager:
    """Prevents text overlap in Manim animations."""
    
    def __init__(self, scene, debug: bool = False):
        self.scene = scene
        self.debug = debug
        self.all_mobjects: List = []
        self.frame_width = config.frame_width
        self.frame_height = config.frame_height
        
    def register_mobject(self, mobject) -> None:
        """Register a mobject for collision tracking."""
        if mobject not in self.all_mobjects:
            self.all_mobjects.append(mobject)
            if self.debug:
                print(f"Registered: {type(mobject).__name__}")
    
    def register_multiple(self, mobjects) -> None:
        """Register multiple mobjects at once."""
        if hasattr(mobjects, '__iter__'):
            for mob in mobjects:
                self.register_mobject(mob)
        else:
            self.register_mobject(mobjects)
    
    def get_bounds(self, mobject) -> Tuple[float, float, float, float]:
        """Get bounding box of mobject (left, right, bottom, top)."""
        try:
            return (
                mobject.get_left()[0], 
                mobject.get_right()[0],
                mobject.get_bottom()[1], 
                mobject.get_top()[1]
            )
        except Exception:
            # Fallback for objects without proper bounds
            center = mobject.get_center()
            return (center[0] - 0.5, center[0] + 0.5, 
                   center[1] - 0.3, center[1] + 0.3)
    
    def check_overlap(self, mob1, mob2, buffer: float = 0.1) -> bool:
        """Check if two mobjects overlap with optional buffer."""
        l1, r1, b1, t1 = self.get_bounds(mob1)
        l2, r2, b2, t2 = self.get_bounds(mob2)
        
        # Add buffer to first mobject's bounds
        l1 -= buffer
        r1 += buffer
        b1 -= buffer
        t1 += buffer
        
        # Check for non-overlap (faster than checking overlap)
        return not (r1 < l2 or l1 > r2 or t1 < b2 or b1 > t2)
    
    def safe_position(self, mobject, preferred_pos: np.ndarray, 
                     strategy: str = "shift") -> np.ndarray:
        """Find safe position avoiding overlaps.
        
        Strategies:
        - 'shift': Move away from overlapping objects
        - 'spiral': Search in expanding spiral pattern
        """
        mobject.move_to(preferred_pos)
        
        # Find overlapping mobjects
        overlaps = [m for m in self.all_mobjects 
                   if m != mobject and self.check_overlap(mobject, m)]
        
        if not overlaps:
            return preferred_pos
        
        if strategy == "shift":
            return self._shift_strategy(mobject, overlaps, preferred_pos)
        else:
            return self._spiral_strategy(mobject, preferred_pos)
    
    def _shift_strategy(self, mobject, overlaps: List, 
                       preferred_pos: np.ndarray) -> np.ndarray:
        """Shift away from overlapping objects."""
        shift = np.zeros(3)
        
        for overlap in overlaps:
            diff = mobject.get_center() - overlap.get_center()
            norm = np.linalg.norm(diff)
            if norm > 0:
                # Normalize and scale by fixed amount
                shift += (diff / norm) * 0.5
        
        new_pos = preferred_pos + shift
        mobject.move_to(new_pos)
        
        # If still overlapping, fall back to spiral
        if any(self.check_overlap(mobject, m) 
               for m in self.all_mobjects if m != mobject):
            return self._spiral_strategy(mobject, preferred_pos)
        
        return new_pos
    
    def _spiral_strategy(self, mobject, center: np.ndarray) -> np.ndarray:
        """Search for safe position in expanding spiral."""
        for radius in np.arange(0.3, 3.0, 0.3):
            for angle in np.arange(0, 2 * np.pi, np.pi / 4):
                test_pos = center + radius * np.array([
                    np.cos(angle), np.sin(angle), 0
                ])
                mobject.move_to(test_pos)
                
                if not any(self.check_overlap(mobject, m) 
                          for m in self.all_mobjects if m != mobject):
                    return test_pos
        
        # Last resort: return shifted position
        return center + np.array([0, 2.0, 0])
'''

# SCENE TEMPLATES AND GUIDANCE


class SceneTemplates:
    """Manages educational video templates and guidance."""
    
    @staticmethod
    def get_progression_guide() -> str:
        """Return the scene progression guide."""
        return """## Educational Video Scene Structure

Create a 2-minute educational video following this progression:

### Scene 1: HOOK (15-20 seconds)
- Introduce what the concept/topic is
- Explain where and why it's used/important
- Create curiosity and relevance

### Scene 2: PHENOMENON (30-40 seconds)
- Show the concept in action with appropriate visuals
- Use animations that best represent the concept
- Ensure scientific/factual accuracy
- Make it visually engaging and clear

### Scene 3: MECHANISM (40-50 seconds)
- Explain how it works in detail
- Break down complex parts into understandable pieces
- Address common misconceptions or confusions
- Use step-by-step visuals when appropriate

### Scene 4: SYNTHESIS (20-30 seconds)
- Summarize the key takeaways (2-3 points)
- End with ONE thought-provoking question or problem
- The question should encourage deeper thinking or application

## Guidelines:
- Choose visuals that best explain the concept
- Use colors meaningfully (consistency helps understanding)
- Transitions should feel natural and FAST (max 0.5s waits)
- Keep narration conversational and clear
- Break long narrations into shorter segments for better sync
- Total duration: approximately 2 minutes

## Voiceover Timing - CRITICAL FOR SYNCHRONIZATION:
- USE MULTIPLE SMALL VOICEOVER BLOCKS instead of one large block
- Each animation should have its own voiceover block for perfect sync
- For Write() animations: Use run_time=tracker.duration (no manual waits needed)
- For other animations: Use run_time=tracker.duration
- AVOID self.wait() inside voiceover blocks - it breaks synchronization
- Keep scene transitions minimal: self.wait(0.5) maximum between scenes

## Correct Pattern:
```python
with self.voiceover(text="First part of explanation") as tracker:
    self.play(Write(title), run_time=tracker.duration)

with self.voiceover(text="Second part continues here") as tracker:
    self.play(Create(circle), run_time=tracker.duration)
```

## WRONG Pattern (causes delays):
```python
with self.voiceover(text="Long explanation covering multiple animations") as tracker:
    self.play(Write(title))
    self.wait(tracker.duration * 0.5)  # This causes delays!
    self.play(Create(circle))
    self.wait(tracker.duration * 0.5)
```

## Collision Avoidance:
IMPORTANT: The CollisionAvoidanceManager class is already provided above. DO NOT create your own.
- Use collision_mgr = CollisionAvoidanceManager(self) at scene start
- Register existing elements: collision_mgr.register_mobject(element)
- Place new elements safely: safe_pos = collision_mgr.safe_position(element, preferred_pos)
- Create new manager for each scene for clean slate
- The provided class uses get_left(), get_right(), etc. NOT get_bounding_box() which is deprecated"""

    @staticmethod
    def get_code_example(class_name: str, voice: str, topic: str) -> str:
        """Generate code structure example with full collision avoidance code."""
        # Get the actual collision avoidance code that will be injected
        collision_code = CollisionAvoidanceCode.get_code()
        
        return f"""```python
import numpy as np
from manim import *
from manim_voiceover import VoiceoverScene
from manim_voiceover.services.openai import OpenAIService

{collision_code}

class {class_name}(VoiceoverScene):
    \"\"\"Educational video explaining {topic}.\"\"\"
    
    def construct(self):
        self.set_speech_service(OpenAIService(voice="{voice}"))
        
        # Scene 1: Hook
        collision_mgr = CollisionAvoidanceManager(self)
        
        title = Text("Title", font_size=48).to_edge(UP)
        collision_mgr.register_mobject(title)
        
        element = Text("Content")
        safe_pos = collision_mgr.safe_position(element, np.array([0, 0, 0]))
        element.move_to(safe_pos)
        collision_mgr.register_mobject(element)
        
        with self.voiceover(text="First part of the narration") as tracker:
            self.play(Write(title), run_time=tracker.duration)
        
        with self.voiceover(text="Second part continues the explanation") as tracker:
            self.play(Write(element), run_time=tracker.duration)
        
        # Minimal transition time
        self.play(FadeOut(*self.mobjects))
        self.wait(0.5)
        
        # Scene 2: Phenomenon (new collision manager)
        collision_mgr = CollisionAvoidanceManager(self)
        # Continue...
```"""

# ============================================================================
# NAME SANITIZATION UTILITIES
# ============================================================================

class NameSanitizer:
    """Handles name sanitization for files and classes."""
    
    @staticmethod
    def sanitize(text: str) -> str:
        """Convert text to valid Python identifier."""
        # Remove special characters, keep alphanumeric and spaces
        sanitized = re.sub(r'[^a-zA-Z0-9\s]', '', text)
        # Replace multiple spaces with single underscore
        sanitized = re.sub(r'\s+', '_', sanitized)
        # Ensure it starts with a letter
        if sanitized and sanitized[0].isdigit():
            sanitized = f'Topic_{sanitized}'
        return sanitized or 'Untitled'
    
    @staticmethod
    def to_class_name(text: str) -> str:
        """Convert to PascalCase class name."""
        sanitized = NameSanitizer.sanitize(text)
        parts = sanitized.split('_')
        return ''.join(word.capitalize() for word in parts if word) + 'Explainer'
    
    @staticmethod
    def to_file_name(text: str) -> str:
        """Convert to snake_case file name."""
        sanitized = NameSanitizer.sanitize(text)
        return sanitized.lower()

# ============================================================================
# VIDEO GENERATOR
# ============================================================================

class VideoGenerator:
    """Handles video generation and execution."""
    
    def __init__(self):
        self.sanitizer = NameSanitizer()
        self.templates = SceneTemplates()
        self.collision_code = CollisionAvoidanceCode()
    
    def create_guidance(self, request: VideoRequest) -> str:
        """Generate guidance for video creation."""
        class_name = self.sanitizer.to_class_name(request.topic)
        
        return f"""LEAP Video Framework Ready

Topic: {request.topic}
Voice: {request.voice.value}
Class Name: {class_name}

{self.templates.get_progression_guide()}

## Your Task:
Create a complete Manim scene that follows the 4-scene structure above.
For the topic "{request.topic}", consider:

1. Hook: What makes {request.topic} interesting or important?
2. Phenomenon: What's the best visual way to show {request.topic} in action?
3. Mechanism: What are the core principles/steps to understand?
4. Synthesis: What's the key insight and what question would make viewers think?

## Code Structure Example:
{self.templates.get_code_example(class_name, request.voice.value, request.topic)}

IMPORTANT: The CollisionAvoidanceManager class shown above is automatically included.
DO NOT write your own collision avoidance implementation.
Use it exactly as shown in the example above."""
    
    def execute(self, request: ExecutionRequest) -> str:
        """Execute the video generation."""
        try:
            # Prepare names and paths
            file_name = self.sanitizer.to_file_name(request.topic)
            class_name = self.sanitizer.to_class_name(request.topic)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Setup directories
            paths = self._setup_directories(file_name, timestamp)
            
            # Process and inject code
            processed_code = self._process_code(
                request.scene_code, class_name, request.topic
            )
            
            # Save script
            self._save_script(processed_code, paths['script'])
            
            # Execute Manim
            return self._execute_manim(
                processed_code, paths, class_name, file_name, timestamp
            )
            
        except Exception as e:
            return self._format_error(e, paths.get('script'))
    
    def _setup_directories(self, file_name: str, timestamp: str) -> Dict[str, Path]:
        """Setup required directories and return paths."""
        scripts_dir = Path.cwd() / Config.SCRIPTS_DIR
        scripts_dir.mkdir(exist_ok=True)
        
        topic_dir = scripts_dir / file_name
        topic_dir.mkdir(exist_ok=True)
        
        output_dir = Path.cwd() / Config.MEDIA_DIR / file_name
        output_dir.mkdir(parents=True, exist_ok=True)
        
        return {
            'script': topic_dir / f"{file_name}_{timestamp}.py",
            'output': output_dir,
            'temp': Path(tempfile.mkdtemp(prefix="leap_mcp_"))
        }
    
    def _process_code(self, code: str, class_name: str, topic: str) -> str:
        """Process and enhance the scene code."""
        # Fix class name
        code = re.sub(
            r'class\s+\w+Explainer\s*\(VoiceoverScene\)',
            f'class {class_name}(VoiceoverScene)',
            code
        )
        
        # Add required imports
        imports = [
            ("from manim import *", "from manim import *\n"),
            ("from manim_voiceover", 
             "from manim_voiceover import VoiceoverScene\n"
             "from manim_voiceover.services.openai import OpenAIService\n"),
            ("import numpy", "import numpy as np\n")
        ]
        
        for check, import_line in imports:
            if check not in code:
                code = import_line + code
        
        # Inject collision avoidance if needed
        if "CollisionAvoidanceManager" not in code:
            code = self._inject_collision_code(code)
        
        # Add header
        header = f'''"""
Generated by LEAP MCP v{Config.VERSION}
Topic: {topic}
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
"""

'''
        return header + code
    
    def _inject_collision_code(self, code: str) -> str:
        """Inject collision avoidance code after imports."""
        collision_code = self.collision_code.get_code()
        
        # Find the last import statement
        import_positions = [
            code.rfind("from manim_voiceover"),
            code.rfind("from manim import"),
            code.rfind("import ")
        ]
        import_end = max(import_positions)
        
        if import_end > -1:
            newline_pos = code.find("\n", import_end)
            if newline_pos > -1:
                return (code[:newline_pos+1] + "\n" + 
                       collision_code + "\n" + 
                       code[newline_pos+1:])
        
        return collision_code + "\n\n" + code
    
    def _save_script(self, code: str, path: Path) -> None:
        """Save the processed script."""
        path.write_text(code, encoding='utf-8')
    
    def _execute_manim(self, code: str, paths: Dict[str, Path], 
                       class_name: str, file_name: str, 
                       timestamp: str) -> str:
        """Execute Manim and return result."""
        temp_script = paths['temp'] / "scene.py"
        temp_script.write_text(code, encoding='utf-8')
        
        output_filename = f"{file_name}_{timestamp}.mp4"
        
        try:
            # Determine quality flag
            quality_map = {
                "low": "-ql",
                "medium": "-qm", 
                "high": "-qh",
                "fourk": "-qk"
            }
            quality_flag = quality_map.get(Config.MANIM_QUALITY, "-qm")
            
            result = subprocess.run(
                [
                    "manim",
                    quality_flag,
                    "-p",  # Preview after rendering
                    "--disable_caching",
                    "-o", output_filename,
                    "--media_dir", str(paths['output']),
                    str(temp_script),
                    class_name
                ],
                capture_output=True,
                text=True,
                timeout=Config.MANIM_TIMEOUT
            )
            
            if result.returncode == 0:
                return self._format_success(
                    paths['output'], output_filename, 
                    paths['script'], class_name
                )
            else:
                return self._format_render_error(
                    result.stderr, paths['script'], class_name
                )
                
        except subprocess.TimeoutExpired:
            return self._format_timeout_error(paths['script'])
        finally:
            # Clean up temp directory
            if paths['temp'].exists():
                shutil.rmtree(paths['temp'], ignore_errors=True)
    
    def _format_success(self, output_dir: Path, output_file: str, 
                       script_path: Path, class_name: str) -> str:
        """Format success message."""
        return f""" Video Created Successfully!

Location: {output_dir}
Video File: {output_file}
Script Saved: {script_path}
Class Name: {class_name}

To re-run this script:
manim -qm -p {script_path} {class_name}"""
    
    def _format_render_error(self, error: str, script_path: Path, 
                            class_name: str) -> str:
        """Format rendering error message."""
        # Extract meaningful error lines
        error_lines = error.split('\n')
        relevant_errors = [line for line in error_lines 
                          if any(keyword in line.lower() 
                                for keyword in ['error', 'exception', 'failed'])]
        error_summary = '\n'.join(relevant_errors[:5]) if relevant_errors else error[:500]
        
        return f""" [Error] Rendering Failed

Error Summary:
{error_summary}

Script saved at: {script_path}

Debug steps:
1. Check the script for syntax errors
2. Verify all Manim objects are imported
3. Ensure class name matches: {class_name}

Manual run command:
manim -qm -p {script_path} {class_name}"""
    
    def _format_timeout_error(self, script_path: Path) -> str:
        """Format timeout error message."""
        return f""" [Error] Video generation timed out (>{Config.MANIM_TIMEOUT}s)

Script saved at: {script_path}

This might be due to:
- Complex animations taking too long
- Infinite loops in the code
- System resource constraints

Try running manually with lower quality:
manim -ql {script_path}"""
    
    def _format_error(self, error: Exception, script_path: Optional[Path]) -> str:
        """Format general error message."""
        error_msg = f""" [Error] Unexpected Error: {str(error)}"""
        
        if script_path and script_path.exists():
            error_msg += f"\n\nScript saved at: {script_path}"
        
        return error_msg

# ============================================================================
# MCP SERVER
# ============================================================================

# Initialize FastMCP server
mcp = FastMCP(Config.SERVER_NAME)

# Initialize video generator
generator = VideoGenerator()

@mcp.tool()
def create_educational_video(topic: str, voice: str = "nova") -> str:
    """
    Create an educational video with structured guidance.
    
    Args:
        topic: The educational topic to explain
        voice: OpenAI TTS voice (alloy, echo, fable, nova, onyx, shimmer)
    
    Returns:
        Guidance and template for creating the video
    """
    try:
        # Validate voice option
        voice_enum = VoiceOption(voice.lower())
        request = VideoRequest(topic=topic, voice=voice_enum)
        return generator.create_guidance(request)
    except ValueError as e:
        return f" Invalid input: {str(e)}"
    except Exception as e:
        return f" Error creating guidance: {str(e)}"

@mcp.tool()
def execute_educational_video(scene_code: str, topic: str = "educational") -> str:
    """
    Execute the educational video code to generate the video.
    
    Args:
        scene_code: Complete Manim scene code
        topic: Topic name for file and class naming
    
    Returns:
        Status message with video location or error details
    """
    try:
        request = ExecutionRequest(scene_code=scene_code, topic=topic)
        return generator.execute(request)
    except ValueError as e:
        return f" Invalid input: {str(e)}"
    except Exception as e:
        return f" Error executing video: {str(e)}"

# Entry point
if __name__ == "__main__":
    mcp.run()
