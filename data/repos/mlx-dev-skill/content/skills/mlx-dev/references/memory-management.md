# MLX Memory Management Reference

## Contents

1. [Memory APIs](#memory-apis)
2. [Critical Gotchas](#critical-gotchas)
3. [Debugging Memory Issues](#debugging-memory-issues)
4. [Best Practices](#best-practices)

---

## Memory APIs

```python
mx.get_active_memory()     # Current usage (bytes)
mx.get_peak_memory()       # Maximum reached (bytes)
mx.get_cache_memory()      # Cached but reusable (bytes)
mx.clear_cache()           # Free cached allocations
mx.set_memory_limit(bytes) # Set ceiling
mx.reset_peak_memory()     # Reset peak counter
```

### Human-Readable Helper

```python
def memory_stats():
    active = mx.get_active_memory() / 1e9
    peak = mx.get_peak_memory() / 1e9
    cache = mx.get_cache_memory() / 1e9
    print(f"Active: {active:.2f}GB, Peak: {peak:.2f}GB, Cache: {cache:.2f}GB")
```

---

## Critical Gotchas

### Memory Is Cached, Not Freed

```python
a = mx.random.normal((10000, 10000))
mx.eval(a)
del a
# Memory still held! Only reused for future allocations

mx.clear_cache()  # Explicitly release
```

### Gradient Accumulation Leaks Without Eval

```python
# WRONG - graph keeps growing
for micro_batch in batches:
    loss, grads = loss_and_grad_fn(model, micro_batch)
    accumulated = tree_map(lambda a, g: a + g, accumulated, grads)

# CORRECT - evaluate each step
for micro_batch in batches:
    loss, grads = loss_and_grad_fn(model, micro_batch)
    accumulated = tree_map(lambda a, g: a + g, accumulated, grads)
    mx.eval(accumulated)
```

### Lazy Loading Only Works With File Paths

```python
# LAZY - arrays loaded on demand
weights = mx.load("model.safetensors")

# NOT LAZY - immediate load
with open("model.safetensors", "rb") as f:
    weights = mx.load(f)  # All arrays evaluated now!
```

---

## Debugging Memory Issues

### Common Symptoms and Causes

| Symptom | Likely Cause |
|---------|--------------|
| OOM during training | Missing mx.eval() in loop |
| Memory grows over epochs | Graph accumulation without eval |
| High cache, low active | Normal caching behavior |
| Can't load model | Model too large for RAM |

### Bundled Memory Utility

Use the bundled `check_memory.py` script:

```bash
# Show current memory stats
uv run python scripts/check_memory.py

# Monitor memory continuously
uv run python scripts/check_memory.py --watch

# Clear cache and show stats
uv run python scripts/check_memory.py --clear
```

### Inline Debugging

```python
import mlx.core as mx

def debug_memory(label=""):
    active = mx.get_active_memory() / 1e9
    peak = mx.get_peak_memory() / 1e9
    cache = mx.get_cache_memory() / 1e9
    print(f"[{label}] Active: {active:.2f}GB, Peak: {peak:.2f}GB, Cache: {cache:.2f}GB")

# Use throughout training
debug_memory("start")
for i, batch in enumerate(dataloader):
    loss = train_step(batch)
    mx.eval(loss, model.parameters())
    if i % 100 == 0:
        debug_memory(f"step {i}")
```

---

## Best Practices

1. **Always mx.eval() at loop boundaries** - Prevents graph accumulation

2. **Use file paths for mx.load()** - Enables lazy loading

3. **Call clear_cache() between major phases** - Reclaims unused memory

4. **Monitor peak memory during development** - Catch leaks early

5. **Use smaller batch sizes if OOM** - MLX can use ~75% of system RAM

6. **Use quantization for large models** - 4-bit reduces memory 4x
