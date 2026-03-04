# MLX Random Number Generation Reference

## Contents

1. [Two Modes](#two-modes)
2. [Critical Gotchas](#critical-gotchas)
3. [Key Splitting Pattern](#key-splitting-pattern)
4. [Best Practices](#best-practices)

---

## Two Modes

### Global Seed (Mutable State)

```python
mx.random.seed(42)
x = mx.random.normal((3, 3))
```

- Simpler for scripts
- State mutates on each call
- Must capture `mx.random.state` in compiled functions

### Explicit Keys (Functional, Reproducible)

```python
key = mx.random.key(42)
x = mx.random.normal((3, 3), key=key)
```

- Deterministic and reproducible
- Requires manual key management
- Better for debugging and testing

---

## Critical Gotchas

### seed() and key() Produce Different Sequences

```python
mx.random.seed(0)
a = mx.random.normal((3,))

key = mx.random.key(0)
b = mx.random.normal((3,), key=key)
# a and b are DIFFERENT!
```

### float16 With Global Seed Can Produce -inf

```python
mx.random.seed(0)
x = mx.random.normal((1, 64, 64, 4), dtype=mx.float16)
x.min()  # Can be -inf!

# Safe: use explicit key
key = mx.random.key(0)
x = mx.random.normal((1, 64, 64, 4), dtype=mx.float16, key=key)
```

**Why this happens:** The Box-Muller transform used for normal distribution sampling can produce extreme values. In float16's limited range (max ~65504), these values overflow to -inf. Explicit keys use a different algorithm that avoids this issue.

### Keys Must Be Split Before Reuse

```python
key = mx.random.key(42)

# WRONG - same result each call
for _ in range(3):
    print(mx.random.uniform(key=key))  # Same value!

# CORRECT - split before each use
for _ in range(3):
    key, subkey = mx.random.split(key)
    print(mx.random.uniform(key=subkey))
```

### Random State in Compiled Functions {#random-state-in-compiled-functions}

Without capturing `mx.random.state`, Dropout produces deterministic (wrong) results:

```python
from functools import partial

state = [model.state, optimizer.state, mx.random.state]  # Include random!

@partial(mx.compile, inputs=state, outputs=state)
def train_step(x, y):
    # Dropout will work correctly
    return loss
```

---

## Key Splitting Pattern

For multiple independent random streams:

```python
key = mx.random.key(42)

# Split for parallel operations
key1, key2, key3 = mx.random.split(key, 3)

# Use each key independently
noise1 = mx.random.normal((100,), key=key1)
noise2 = mx.random.normal((100,), key=key2)
indices = mx.random.randint(0, 10, (50,), key=key3)
```

---

## Best Practices

1. **Use explicit keys for reproducibility** - Easier to debug

2. **Split keys, don't reuse** - Same key = same values

3. **Capture random state in compiled functions** - For correct Dropout

4. **Avoid float16 with global seed** - Use explicit key instead

5. **For data augmentation, split per sample**:
   ```python
   keys = mx.random.split(base_key, batch_size)
   augmented = mx.vmap(augment)(data, keys)
   ```
