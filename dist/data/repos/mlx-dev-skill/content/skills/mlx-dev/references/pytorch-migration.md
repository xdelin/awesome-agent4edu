# PyTorch to MLX Migration Reference

## Contents

1. [Key Differences](#key-differences)
2. [Weight Conversion](#weight-conversion)
3. [API Mapping](#api-mapping)
4. [Serialization](#serialization)

---

## Key Differences

| Aspect | PyTorch | MLX |
|--------|---------|-----|
| Device management | `tensor.to('cuda')` | Not needed (unified memory) |
| Conv format | NCHW | NHWC |
| Forward method | `forward()` | `__call__()` |
| Eager/Lazy | Eager by default | Lazy by default |
| Gradient sync | Automatic | Manual `mx.eval()` |

---

## Weight Conversion

### Basic Conversion

```python
import torch
import mlx.core as mx

# Load PyTorch weights
pt_weights = torch.load("model.pt", map_location="cpu")

# Convert to MLX
mlx_weights = {k: mx.array(v.numpy()) for k, v in pt_weights.items()}
```

### Conv2d Weight Format Change

PyTorch: (out_channels, in_channels, kH, kW)
MLX: (out_channels, kH, kW, in_channels)

```python
def convert_conv_weights(pt_weights):
    mlx_weights = {}
    for k, v in pt_weights.items():
        if "conv" in k and "weight" in k and len(v.shape) == 4:
            # OIHW -> OHWI
            mlx_weights[k] = mx.array(v.numpy().transpose(0, 2, 3, 1))
        else:
            mlx_weights[k] = mx.array(v.numpy())
    return mlx_weights
```

### BatchNorm Conversion

```python
# PyTorch names → MLX names
name_mapping = {
    "running_mean": "running_mean",
    "running_var": "running_var",
    "weight": "weight",  # gamma
    "bias": "bias",      # beta
}
```

---

## API Mapping

### Core Operations

| PyTorch | MLX |
|---------|-----|
| `torch.tensor()` | `mx.array()` |
| `tensor.cuda()` | Not needed |
| `tensor.cpu()` | Not needed |
| `torch.no_grad()` | Not needed (lazy) |
| `tensor.detach()` | `mx.stop_gradient()` |
| `tensor.item()` | `arr.item()` |

### Neural Network Layers

| PyTorch | MLX |
|---------|-----|
| `nn.Linear` | `nn.Linear` |
| `nn.Conv2d` | `nn.Conv2d` (NHWC!) |
| `nn.LayerNorm` | `nn.LayerNorm` |
| `nn.BatchNorm2d` | `nn.BatchNorm` |
| `nn.Dropout` | `nn.Dropout` |
| `nn.MultiheadAttention` | `nn.MultiHeadAttention` |

### Tensor Operations

| PyTorch | MLX |
|---------|-----|
| `torch.gather()` | `mx.take_along_axis()` |
| `torch.scatter_add_()` | `arr.at[idx].add()` |
| `torch.einsum()` | `mx.einsum()` |
| `tensor.view()` | `arr.reshape()` |
| `tensor.permute()` | `mx.transpose()` |
| `tensor.contiguous()` | Not needed |

### Training Operations

| PyTorch | MLX |
|---------|-----|
| `loss.backward()` | `nn.value_and_grad()` |
| `optimizer.step()` | `optimizer.update()` |
| `optimizer.zero_grad()` | Not needed |
| `model.train()` | Pass training flag |
| `model.eval()` | Pass training flag |

---

## Serialization

### MLX Safetensors Have Incompatible Metadata

Files saved by MLX have `format: "mlx"` instead of `"pt"`, which can prevent HuggingFace transformers from loading them.

### HuggingFace Compatibility Workaround

To create weights compatible with HuggingFace:

```python
import safetensors.torch as st_torch

# Convert MLX weights to numpy, then save with safetensors.torch
weights_np = {k: np.array(v.astype(mx.float32)) for k, v in model.parameters().items()}
weights_pt = {k: torch.from_numpy(v) for k, v in weights_np.items()}
st_torch.save_file(weights_pt, "model_hf_compatible.safetensors")
```

Alternatively, strip the metadata when loading in HuggingFace by using `safetensors` directly instead of the transformers auto-loader.

### Saving Weights

```python
# Single array
mx.save("array.npy", arr)

# Multiple arrays (safetensors format)
mx.savez("weights.safetensors", **{"layer1": w1, "layer2": w2})

# Model weights
model.save_weights("model.safetensors")
```

### Loading Weights

```python
# Single array
arr = mx.load("array.npy")

# Multiple arrays
weights = mx.load("weights.safetensors")

# Model weights
model.load_weights("model.safetensors")
```
