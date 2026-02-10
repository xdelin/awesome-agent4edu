# mlx-dev-skill

A Claude Code skill for writing correct, idiomatic [Apple MLX](https://github.com/ml-explore/mlx) code on Apple Silicon.

## Installation

### Option 1: Plugin Marketplace (Recommended)

In Claude Code, run:

```
/plugin marketplace add tkwn2080/mlx-dev-skill
```

The skill will be automatically available after adding the marketplace.

### Option 2: Manual Installation

Clone and copy to your personal skills directory:

```bash
git clone https://github.com/tkwn2080/mlx-dev-skill.git
cp -r mlx-dev-skill/skills/mlx-dev ~/.claude/skills/
```

Then restart Claude Code.

## What This Skill Provides

When working with MLX, Claude will automatically:

- Apply correct lazy evaluation patterns with `mx.eval()` at loop boundaries
- Use proper array indexing (lists must be `mx.array`, slices create copies)
- Use NHWC format for Conv2d (not NCHW like PyTorch)
- Override `__call__()` not `forward()` for nn.Module
- Avoid float64 on GPU (CPU-only)
- Capture all mutable state in `mx.compile()`
- Apply proper memory management patterns

## Coverage

- **Array Operations**: Indexing, slicing, `at[]` syntax, gather/scatter
- **Neural Networks**: Layer equivalents, weight formats, quantization
- **Compilation**: `mx.compile()` patterns, state capture, shapeless mode
- **Memory**: Debugging tools, cache management, leak detection
- **Data Types**: GPU support table, bfloat16 handling, conversion patterns
- **Random**: Seed vs key, splitting patterns, compiled function state
- **Gradients**: `value_and_grad`, custom vjp, control flow
- **PyTorch Migration**: Weight conversion, API mapping, format changes
- **Error Decoding**: Common errors mapped to solutions

## Structure

```
mlx-dev-skill/
├── .claude-plugin/
│   └── plugin.json             # Plugin metadata
├── skills/
│   └── mlx-dev/
│       ├── SKILL.md            # Main skill entry point
│       ├── references/
│       │   ├── array-indexing.md
│       │   ├── compilation.md
│       │   ├── dtypes.md
│       │   ├── error-decoder.md
│       │   ├── gradients.md
│       │   ├── memory-management.md
│       │   ├── neural-networks.md
│       │   ├── pytorch-migration.md
│       │   └── random.md
│       └── scripts/
│           └── check_memory.py
├── README.md
└── LICENSE
```

## Memory Debugging Utility

The skill includes a memory debugging script:

```bash
# Show current memory stats
uv run python ~/.claude/skills/mlx-dev/scripts/check_memory.py

# Monitor continuously
uv run python ~/.claude/skills/mlx-dev/scripts/check_memory.py --watch

# Log to CSV for analysis
uv run python ~/.claude/skills/mlx-dev/scripts/check_memory.py --watch --log memory.csv
```

## Requirements

- Apple Silicon Mac (M1/M2/M3/M4)
- MLX installed (`uv add mlx`)
- Claude Code CLI

## License

MIT
