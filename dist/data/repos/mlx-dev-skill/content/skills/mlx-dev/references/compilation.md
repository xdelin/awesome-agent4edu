# MLX Compilation Reference

## Contents

1. [Basic Usage](#basic-usage)
2. [Critical Gotchas](#critical-gotchas)
3. [Mutable State Capture](#mutable-state-capture)
4. [Recompilation Triggers](#recompilation-triggers)
5. [Decision Tree](#decision-tree)

---

## Basic Usage

```python
@mx.compile
def gelu(x):
    return x * (1 + mx.erf(x / math.sqrt(2))) / 2
# ~5x speedup on M1 Max
```

---

## Critical Gotchas

### Print Statements Crash Compiled Functions

```python
@mx.compile
def fn(x):
    print(x)  # CRASH - placeholder has no data during tracing
    return mx.exp(x)

# Debug with:
mx.disable_compile()
fn(mx.array(5.0))
```

### Random State Must Be Captured

See [random.md#random-state-in-compiled-functions](random.md#random-state-in-compiled-functions) for complete documentation on capturing random state in compiled functions.

### String Decoding Triggers Recompilation

```python
# WRONG - recompiles every iteration, eventual OOM!
for step in range(1000):
    text_decoded = some_bytes.decode("utf-8")
    loss = compiled_step(x, text_decoded, y)

# CORRECT - use literal or pre-decoded string
text = some_bytes.decode("utf-8")
for step in range(1000):
    loss = compiled_step(x, text, y)
```

### Shape-Dependent Code Breaks with shapeless=True

```python
# WRONG with shapeless=True
def fn(x):
    return x.reshape(x.shape[0] * x.shape[1], -1)

# CORRECT - shape-agnostic
def fn(x):
    return x.flatten(0, 1)
```

---

## Using shapeless=True

Use `shapeless=True` when:
- Input shapes vary between calls (variable batch size, sequence length)
- You want to avoid recompilation on shape changes

```python
@partial(mx.compile, shapeless=True)
def flexible_fn(x):
    return x.flatten(0, 1)  # Works with any batch size
```

**Requirements:**
- Avoid shape-dependent code (no `x.shape[0] * x.shape[1]`)
- Use shape-agnostic operations (`.flatten()` instead of `.reshape(N, M)`)

**Trade-off:** Slightly less optimized kernels but no recompilation overhead.

---

## Mutable State Capture

For functions with mutable state (model parameters, optimizer, random state):

```python
from functools import partial

state = [model.state, optimizer.state, mx.random.state]

@partial(mx.compile, inputs=state, outputs=state)
def train_step(x, y):
    loss, grads = nn.value_and_grad(model, loss_fn)(model, x, y)
    optimizer.update(model, grads)
    return loss
```

**What must be captured:**
- `model.state` - Model parameters that get updated
- `optimizer.state` - Optimizer momentum/velocity buffers
- `mx.random.state` - For Dropout, random augmentation, etc.

---

## Recompilation Triggers

Compiled functions recompile when:

| Trigger | Example |
|---------|---------|
| Input shapes change | Batch size varies |
| Input dtypes change | float32 → float16 |
| Number of inputs changes | Optional arguments |
| New Python objects created inside | String decoding in loop |

**Symptoms of excessive recompilation:**
- Slower than uncompiled version
- Growing memory usage (OOM)
- Inconsistent timing

---

## Decision Tree

```
Is function called more than once?
  NO → Don't compile
  YES ↓

Does function have mutable state (model, optimizer, random)?
  YES → Use @partial(mx.compile, inputs=state, outputs=state)
  NO ↓

Does function have data-dependent control flow?
  YES → Restructure with mx.where or don't compile
  NO ↓

Does function need to work with varying shapes?
  YES → Consider shapeless=True, avoid shape-dependent code
  NO → Simple @mx.compile is fine
```
