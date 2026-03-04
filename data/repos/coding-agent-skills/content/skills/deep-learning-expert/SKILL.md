---
name: deep-learning-expert
description: Use this skill when user needs to design deep learning models, perform training optimization, interpret papers, or solve engineering practice issues. Trigger keywords: deep learning, neural networks, model training, overfitting, gradient vanishing, PyTorch, TensorFlow, pretrained models, fine-tuning, model deployment, inference optimization, CV (Computer Vision), NLP (Natural Language Processing), reinforcement learning. Applicable to full-cycle scenarios of deep learning modeling, training, debugging, and deployment.
---

# Deep Learning Modeling and Optimization

## Description
Provide deep learning model design, training optimization, problem diagnosis, and engineering practice guidance to ensure models meet requirements for accuracy, efficiency, and deployability.

## When to Use
- User requests "design a deep learning model" or "how to train XX task"
- User asks about training issues (overfitting, underfitting, gradient vanishing/explosion, NaN, slow convergence)
- User needs paper reproduction or cutting-edge technology analysis (Transformer, Diffusion, ViT, etc.)
- User seeks framework usage recommendations (PyTorch, TensorFlow, JAX) or code implementation
- User asks about best practices in CV, NLP, reinforcement learning domains
- User needs model compression, quantization, deployment, or inference optimization solutions
- User requests data processing, data augmentation, or sampling strategy recommendations
- User compares different model architectures or training strategies

## When NOT to Use
- User only needs traditional machine learning methods (e.g., decision trees, SVM, linear regression) without deep learning
- User's problem is pure data analysis or statistics without model training
- User only needs data preprocessing or feature engineering without model design
- User's problem is algorithm competition or pure programming without deep learning modeling
- User needs MLOps platform setup or team management rather than specific technical solutions

## Input
```typescript
{
  task: {
    type: string                  // Task type (classification/detection/segmentation/generation/sequence-modeling/RL)
    domain: string                // Application domain (CV/NLP/audio/multimodal/time-series)
    specificTask?: string         // Specific task (e.g., image-classification → CIFAR-10)
  }
  data: {
    size?: string                 // Dataset size (e.g., "10k samples")
    quality?: string              // Data quality description (noise, annotation quality)
    distribution?: string         // Class distribution (balanced/long-tail)
    labelAvailability?: string    // Labeling status (fully-supervised/weakly-supervised/unsupervised)
    augmentationNeeded?: boolean  // Whether data augmentation is needed
  }
  constraints: {
    computeResources?: string     // Compute resources (single-GPU/multi-GPU/TPU/CPU)
    latencyRequirement?: string   // Latency requirement (e.g., "inference <50ms")
    accuracyTarget?: string       // Accuracy target (e.g., "top-1 acc >90%")
    modelSizeLimit?: string       // Model size limit (e.g., "<100MB")
    budget?: string               // Budget constraints (affects pretrained model selection/compute resources)
  }
  challenges?: string[]           // Known challenges (few-shot/domain-shift/real-time/long-tail)
  existingSetup?: {
    framework?: string            // Existing framework (PyTorch/TensorFlow/JAX)
    baseModel?: string            // Existing base model
    trainingIssues?: string[]     // Current training issues
  }
}
```

## Output
```typescript
{
  modelDesign: {
    architecture: string          // Recommended architecture (CNN/Transformer/RNN/GNN/Diffusion/Hybrid)
    rationale: string             // Architecture selection rationale (inductive bias/parameter count/computational complexity)
    pretrainedModel?: {
      name: string                // Pretrained model (e.g., ResNet-50, BERT-base, ViT-B/16)
      source: string              // Source (HuggingFace/TorchVision/OpenAI)
      adaptation: string          // Adaptation method (full fine-tuning/LoRA/Adapter/freeze partial layers)
    }
    tradeoff: string              // Train from scratch vs fine-tuning vs transfer learning tradeoff analysis
    codeFramework?: string        // PyTorch/TensorFlow code framework (if applicable)
  }
  trainingStrategy: {
    dataProcessing: {
      normalization: string       // Normalization method
      augmentation?: string[]     // Data augmentation strategies
      samplingStrategy?: string   // Sampling strategy (balanced sampling/hard negative mining)
      batchSize: string           // Batch size recommendation (considering memory and convergence)
    }
    optimizer: {
      type: string                // Optimizer type (SGD/AdamW/Lion)
      learningRate: string        // Initial learning rate
      scheduler?: string          // Learning rate scheduler (Cosine/StepLR/ReduceLROnPlateau)
      weightDecay?: string        // Weight decay
    }
    trainingTechniques: {
      gradientClipping?: string   // Gradient clipping threshold
      mixedPrecision?: boolean    // Mixed precision training (FP16/BF16)
      gradientAccumulation?: number // Gradient accumulation steps
      ema?: boolean               // Exponential moving average
    }
    regularization?: string[]     // Regularization methods (Dropout/Batch Norm/Label Smoothing)
  }
  debuggingGuidance: {
    visualization: string[]       // Visualization recommendations (loss curves/weight distributions/activation heatmaps)
    checkpoints: string[]         // Checkpoints (weight norms/gradient norms/data statistics)
    commonIssues: {
      issue: string               // Issue type (overfitting/underfitting/NaN/non-convergence)
      diagnosis: string           // Diagnosis method
      solution: string[]          // Solutions
    }[]
  }
  engineeringPractice: {
    modelCompression?: {
      pruning?: string            // Pruning strategy
      quantization?: string       // Quantization method (INT8/FP16)
      distillation?: string       // Knowledge distillation
    }
    inferenceOptimization?: {
      format: string[]            // Inference formats (ONNX/TensorRT/CoreML)
      batching?: string           // Batching strategy
      caching?: string            // Caching strategy
    }
    distributedTraining?: {
      strategy: string            // Distributed strategy (DDP/DeepSpeed/FSDP)
      configuration: string       // Configuration recommendations
    }
  }
  expectedPerformance?: {
    accuracy?: string             // Expected accuracy
    trainingTime?: string         // Expected training time
    inferenceLatency?: string     // Expected inference latency
  }
}
```

## Execution Workflow

Copy the following checklist before starting, and explicitly mark status after completing each step.

### Step 1: Task Requirement Analysis
- Clarify task type (classification/detection/segmentation/generation/sequence-modeling/RL)
- Confirm application domain (CV/NLP/audio/multimodal/time-series)
- Quantify data situation (size, quality, labeling, distribution)
- Identify constraints (compute resources, latency, accuracy, model size)
- Extract core challenges (few-shot, domain shift, long-tail distribution, real-time)

**Feedback Loop**: If critical information is missing (e.g., data size, compute resources), MUST ask user to provide details. Avoid designing based on assumptions.

### Step 2: Model Architecture Selection
- List candidate architectures (CNN/Transformer/RNN/GNN/Diffusion/hybrid)
- Explain each architecture's inductive bias (e.g., CNN's locality, Transformer's long-range dependencies)
- Compare parameter count, computational complexity, data requirements
- Assess pretrained model availability (HuggingFace/TorchVision/OpenAI)
- Tradeoff train from scratch vs fine-tuning vs transfer learning

**Architecture Selection Principles**:
- **Small data (<10k)**: Prioritize transfer learning or fine-tuning pretrained models
- **Medium data (10k-100k)**: Fine-tuning + data augmentation
- **Large data (>100k)**: Can consider training from scratch or light fine-tuning
- **Real-time requirements**: Prioritize lightweight architectures (MobileNet/EfficientNet/DistilBERT)

**Feedback Loop**: If user is indecisive about multiple architectures, provide detailed comparison of 2-3 options (accuracy/speed/resource consumption).

### Step 3: Training Strategy Design
- Data processing: Normalization methods (ImageNet stats/adaptive), augmentation strategies (AutoAugment/RandAugment)
- Optimizer selection: SGD (vision tasks) vs AdamW (language/multimodal tasks)
- Learning rate scheduling: Warmup + Cosine Annealing (recommended) / StepLR / ReduceLROnPlateau
- Training techniques: Gradient clipping (prevent explosion), mixed precision (accelerate training), gradient accumulation (simulate large batch)
- Regularization: Dropout, weight decay, label smoothing

**Training Recipe Template**:
```python
# Example: Image classification training recipe
optimizer = AdamW(model.parameters(), lr=1e-3, weight_decay=0.01)
scheduler = CosineAnnealingLR(optimizer, T_max=epochs)
scaler = GradScaler()  # Mixed precision
# Data augmentation: RandomCrop + RandomHorizontalFlip + ColorJitter
```

**Feedback Loop**: If user reports training instability (loss oscillation/NaN), immediately check learning rate, data normalization, gradient norms.

### Step 4: Experiment Execution and Monitoring
- Establish baseline experiment (simplest viable model)
- Set monitoring metrics (train/val loss, accuracy, gradient norm, learning rate)
- Configure visualization tools (TensorBoard/WandB)
- Save key checkpoints (best model, last model, every N epochs)

**Monitoring Checklist**:
- [ ] Training loss smoothly decreases (no oscillation)
- [ ] Val loss and train loss gap is reasonable (<2x)
- [ ] Gradient norm is stable (no explosion/vanishing)
- [ ] Learning rate decays as expected
- [ ] Data loading has no bottleneck (GPU utilization >80%)

**Feedback Loop**: If val loss doesn't decrease early, check data pipeline (label errors/normalization issues/excessive augmentation).

### Step 5: Problem Diagnosis and Tuning
Diagnose and resolve common issues based on training performance:

#### Overfitting (high train acc, low val acc)
- **Diagnosis**: Train/val loss curves diverge, val accuracy stops improving
- **Solutions**:
  - Enhance data augmentation (AutoAugment/MixUp/CutMix)
  - Increase Dropout rate or weight decay
  - Use Early Stopping
  - Reduce model capacity (fewer layers/channels)
  - Acquire more training data

#### Underfitting (both train/val acc low)
- **Diagnosis**: Training loss doesn't decrease for long, accuracy below random guess
- **Solutions**:
  - Increase model capacity (deeper/wider network)
  - Reduce regularization strength (lower Dropout/weight decay)
  - Check data quality (annotation errors/class imbalance)
  - Increase training epochs or adjust learning rate

#### NaN/Inf (training interruption)
- **Diagnosis**: Loss suddenly becomes NaN, gradient explosion
- **Solutions**:
  - Reduce learning rate (10x smaller)
  - Enable gradient clipping (clip_grad_norm=1.0)
  - Check input data (contains NaN/Inf)
  - Add numerical stability operations (e.g., add epsilon in log)

#### Slow Convergence
- **Diagnosis**: Loss decreases slowly, training time too long
- **Solutions**:
  - Increase learning rate (careful, may be unstable)
  - Use more aggressive learning rate schedule (e.g., OneCycleLR)
  - Enable mixed precision training (2-3x speedup)
  - Optimize data loading (multi-process, prefetch)

**Feedback Loop**: If multiple tuning attempts still don't improve, revisit Step 2 (architecture selection) or Step 1 (task definition) for reasonableness.

### Step 6: Model Compression and Deployment Optimization (if applicable)
- Assess compression needs (model size/inference latency meets constraints)
- Choose compression strategy:
  - **Pruning**: Structured pruning (reduce channels/layers) vs unstructured pruning (reduce parameters)
  - **Quantization**: PTQ (post-training quantization) vs QAT (quantization-aware training), INT8/FP16
  - **Distillation**: Use large model (teacher) to train small model (student)
- Convert inference format (ONNX/TensorRT/CoreML)
- Optimize inference pipeline (batching/operator fusion/memory optimization)

**Compression Tradeoffs**:
- **Pruning**: Medium accuracy drop (1-3%), high compression ratio (2-5x)
- **Quantization**: Small accuracy drop (<1%), medium compression ratio (2-4x)
- **Distillation**: Controllable accuracy drop, compression ratio depends on student model size

**Feedback Loop**: If accuracy drop after compression exceeds expectation (>5%), consider relaxing compression ratio or using distillation compensation.

### Step 7: Performance Validation and Documentation
- Validate final performance on test set (accuracy, F1, mAP, etc.)
- Measure inference latency (single sample/batching)
- Record model version, training configuration, dataset version
- Generate Model Card: purpose, performance, limitations, ethical considerations

**Performance Benchmark Comparison**:
- Compare with public benchmarks (SOTA)
- Compare with simple baseline (e.g., linear model)
- Analyze failure cases (which classes/samples are error-prone)

**Feedback Loop**: If performance doesn't meet target and no obvious improvement space, discuss with user whether to adjust targets or increase resource investment.

## Failure Handling

### Insufficient or Poor Quality Data
- **Symptom**: User provides extremely small data (<1000 samples) or poor annotation quality
- **Action**: 
  - Recommend data augmentation + pretrained model fine-tuning
  - Suggest active learning or semi-supervised learning
  - If budget allows, suggest outsourcing annotation or synthetic data generation

### Insufficient Compute Resources
- **Symptom**: User only has CPU or single low-end GPU, but task requires large model
- **Action**:
  - Recommend lightweight architectures (MobileNet/DistilBERT/TinyLlama)
  - Suggest using pretrained models + freezing partial layers
  - Provide cloud platform solutions (Colab/Kaggle/AWS SageMaker)

### Training Persistently Fails to Converge
- **Symptom**: Trying multiple learning rates/optimizers still doesn't converge
- **Action**:
  - Check data pipeline (visualize batch samples)
  - Simplify model (start from minimal viable model)
  - Check labels are correct (cross-validate annotation quality)

### Post-Deployment Performance Falls Short
- **Symptom**: Model performs well on test set but poorly in production
- **Action**:
  - Analyze train/production data distribution difference (domain shift)
  - Use domain adaptation techniques (DANN/CORAL)
  - Collect production data for incremental training

### Over-Reliance on Pretrained Models
- **Symptom**: User unconditionally uses large pretrained models, ignoring task characteristics
- **Action**: Explain pretrained model applicable scenarios, compare cost/benefit of training from scratch, recommend moderate solutions