# MLX Neural Networks Reference

## Contents

1. [Critical Differences from PyTorch](#critical-differences-from-pytorch)
2. [Layer Equivalents](#layer-equivalents)
3. [Model Definition Pattern](#model-definition-pattern)
4. [Quantization](#quantization)

---

## Critical Differences from PyTorch

### Conv2d Uses NHWC Format

PyTorch uses NCHW, MLX uses NHWC:

```python
# PyTorch: (N, C, H, W)
x_torch = torch.randn(1, 3, 224, 224)

# MLX: (N, H, W, C)
x_mlx = mx.array(x_torch.numpy().transpose(0, 2, 3, 1))

# Weight conversion:
# PyTorch: (out, in, kH, kW)
# MLX: (out, kH, kW, in)
```

### Override __call__(), Not forward()

```python
class MyModel(nn.Module):
    def __call__(self, x):  # Not forward()!
        return self.layer(x)
```

### No dtype Parameter in Layer Constructors

```python
# PyTorch style - NOT supported:
# layer = nn.Linear(10, 10, dtype=torch.bfloat16)

# MLX pattern:
layer = nn.Linear(10, 10)
layer.set_dtype(mx.bfloat16)
```

### BatchNorm Running Stats

In versions before v0.30, BatchNorm's running_mean and running_var didn't update inside `value_and_grad()`. Fixed in v0.30+.

---

## Layer Equivalents

| PyTorch | MLX |
|---------|-----|
| `nn.Linear` | `nn.Linear` |
| `nn.Conv2d` | `nn.Conv2d` (NHWC!) |
| `nn.LayerNorm` | `nn.LayerNorm` |
| `nn.BatchNorm2d` | `nn.BatchNorm` |
| `nn.Dropout` | `nn.Dropout` |
| `nn.MultiheadAttention` | `nn.MultiHeadAttention` |
| `nn.LSTM` | `nn.LSTM` |
| `nn.GRU` | `nn.GRU` |
| `nn.Embedding` | `nn.Embedding` |

---

## Model Definition Pattern

```python
import mlx.core as mx
import mlx.nn as nn

class TransformerBlock(nn.Module):
    def __init__(self, dims: int, num_heads: int):
        super().__init__()
        self.attention = nn.MultiHeadAttention(dims, num_heads)
        self.norm1 = nn.LayerNorm(dims)
        self.norm2 = nn.LayerNorm(dims)
        self.ffn = nn.Sequential(
            nn.Linear(dims, dims * 4),
            nn.GELU(),
            nn.Linear(dims * 4, dims)
        )

    def __call__(self, x, mask=None):
        # Pre-norm architecture
        h = self.norm1(x)
        h = self.attention(h, h, h, mask=mask)
        x = x + h

        h = self.norm2(x)
        h = self.ffn(h)
        return x + h
```

---

## Quantization

```python
nn.quantize(model, group_size=64, bits=4)
# ValueError if input_dims not divisible by 64
```

**Requirements:**
- Input dimensions must be divisible by group_size
- Common group sizes: 32, 64, 128
- Common bit widths: 4, 8

### Quantization Trade-offs

| Bits | Memory Reduction | Quality | Use Case |
|------|------------------|---------|----------|
| 4-bit | 4x smaller | Good for inference | LLM deployment |
| 8-bit | 2x smaller | Near-lossless | General models |

**group_size selection:**
- 32: More accurate, slightly larger file size
- 64: Balanced (default, recommended)
- 128: Smaller files, potential quality loss

**Saving/Loading Quantized Models:**
```python
# Save
model.save_weights("model_q4.safetensors")

# Load (weights stay quantized)
model.load_weights("model_q4.safetensors")
```
