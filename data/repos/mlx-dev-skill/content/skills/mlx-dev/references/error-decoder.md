# MLX Error Messages Decoded

## Contents

1. [Common Errors](#common-errors)
2. [Silent Failures](#silent-failures)
3. [Debugging Tips](#debugging-tips)

---

## Common Errors

| Error Message | Actual Cause | Solution |
|---------------|--------------|----------|
| `Cannot index mlx array using the given type` | Using Python list for indexing | Wrap with `mx.array([...])`. See [array-indexing.md](array-indexing.md) |
| `Slice indices must be integers or None` | MLX array as slice bound | Use `.item()` to get Python int |
| `Primitive's vjp not implemented` | Non-differentiable operation | Restructure (use concatenate, where). See [gradients.md](gradients.md) |
| `ValueError: Gradients not implemented` | Scatter operations in gradient path | Use concatenate/where instead of indexed assignment |
| `ValueError: must be non-empty array` | Wrong shapes in vmap | Check input dimensions match vmap axes |
| `mach-o incompatible architecture` | Using x86 Python via Rosetta | Install ARM Python |
| `kIOGPUCommandBufferCallbackErrorOutOfMemory` | Model too large or memory leak | Reduce batch, check eval placement. See [memory-management.md](memory-management.md) |
| `Unable to allocate X bytes` | Memory fragmentation | Call `clear_cache()`. See [memory-management.md](memory-management.md) |
| `Unable to load kernel` | Missing Metal kernel | Update MLX version |
| `RuntimeError` on float64 GPU op | float64 not supported on GPU | Use `stream=mx.cpu` or cast to float32. See [dtypes.md](dtypes.md) |
| Metal kernel timeout | Very long single operation | Break into smaller chunks |
| `format: "mlx"` load failure | MLX safetensors vs HuggingFace | Re-export without MLX metadata. See [pytorch-migration.md](pytorch-migration.md) |

---

## Silent Failures

These produce no error but give wrong results:

### Out-of-Bounds Indexing

```python
a = mx.ones([4, 4])
a[0, 100]  # Returns garbage - no IndexError!
```

**Solution:** Validate indices manually before indexing.

### Duplicate Index Updates

```python
a = mx.array([1, 2, 3])
a[[0, 0]] = mx.array([4, 5])  # a[0] could be 4 OR 5!
```

**Solution:** Use `at[].add()` for accumulation, avoid duplicate indices for assignment.

### Dropout in Compiled Functions Without Random State

```python
@mx.compile  # Missing state capture!
def forward(x):
    return nn.Dropout(0.5)(x)  # Returns same pattern every call
```

**Solution:** Capture `mx.random.state` in compiled functions. See [random.md#random-state-in-compiled-functions](random.md#random-state-in-compiled-functions) for complete documentation.

### Integer Overflow in Statistics

```python
a = mx.array([0, 100, 200, 123], dtype=mx.uint8)
mx.mean(a)  # Returns 41.75 - WRONG!
```

**Solution:** Cast to float32 before computing statistics.

### bfloat16 Misinterpretation

```python
from ml_dtypes import bfloat16
x = np.array(1., dtype=bfloat16)
mx.array(x)  # Returns complex64!
```

**Solution:** Convert through float32: `mx.array(x.astype(np.float32), dtype=mx.bfloat16)`

---

## Debugging Tips

### Disable Compilation for Debugging

```python
mx.disable_compile()
# Now print() works in formerly-compiled functions
fn(mx.array(5.0))
mx.enable_compile()
```

### Memory Debugging

```python
def debug_memory(label=""):
    active = mx.get_active_memory() / 1e9
    peak = mx.get_peak_memory() / 1e9
    print(f"[{label}] Active: {active:.2f}GB, Peak: {peak:.2f}GB")

debug_memory("before")
# ... operation ...
debug_memory("after")
```

### Force Evaluation for Debugging

```python
# Force eval to see actual values
result = some_operation(x)
mx.eval(result)
print(result)  # Now shows computed values
```

### Check Shapes Throughout Pipeline

```python
def debug_shape(x, label=""):
    print(f"[{label}] shape: {x.shape}, dtype: {x.dtype}")
    return x
```
