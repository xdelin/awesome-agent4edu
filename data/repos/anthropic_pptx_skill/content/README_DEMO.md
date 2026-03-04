# Anthropic Skills Demo for PPTX

This directory contains the setup for using the `anthropics/skills` repository's PPTX capabilities.

## Status
- **Analysis & Extraction**: ✅ Operational (using `markitdown` and `ooxml/scripts/unpack.py`)
- **Creation (HTML to PPTX)**: ⚠️ Partially blocked (requires `playwright` browser installation which failed due to network).
- **Creation (Python)**: ✅ Verified via `make_sample_ppt.py`.

## Files
- `sample.pptx`: A sample presentation created locally.
- `unpacked_sample/`: The raw XML structure of the sample presentation (result of `unpack.py`).

## How to use
1. **Analyze a PPT**:
   ```bash
   python -m markitdown sample.pptx
   ```

2. **Unpack a PPT for editing**:
   ```bash
   python ooxml/scripts/unpack.py sample.pptx output_dir
   ```

3. **Repack after editing** (if needed):
   ```bash
   python ooxml/scripts/pack.py output_dir new_presentation.pptx
   ```
