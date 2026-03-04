---
name: deep-learning-expert
description: 当用户需要设计深度学习模型、进行训练调优、解读论文或解决工程实践问题时使用此技能。触发关键词：深度学习、神经网络、模型训练、过拟合、梯度消失、PyTorch、TensorFlow、预训练模型、微调、模型部署、推理优化、CV(计算机视觉)、NLP(自然语言处理)、强化学习。适用于深度学习建模、训练、调试和部署的全流程场景。
---

# 深度学习建模与优化

## 描述
提供深度学习模型设计、训练调优、问题诊断和工程实践指导，确保模型在准确性、效率和可部署性上满足需求。

## 何时使用
- 用户请求"设计一个深度学习模型"或"如何训练XX任务"
- 用户询问训练问题(过拟合、欠拟合、梯度消失/爆炸、NaN、收敛慢)
- 用户需要论文复现或前沿技术解析(Transformer、Diffusion、ViT等)
- 用户寻求框架使用建议(PyTorch、TensorFlow、JAX)或代码实现
- 用户询问 CV、NLP、强化学习等领域的最佳实践
- 用户需要模型压缩、量化、部署或推理优化方案
- 用户请求数据处理、数据增强或采样策略建议
- 用户对比不同模型架构或训练策略的优劣

## 何时不使用
- 用户只需要传统机器学习方法(如决策树、SVM、线性回归)，不涉及深度学习
- 用户的问题是纯数据分析或统计学，不涉及模型训练
- 用户只需要数据预处理或特征工程，不涉及模型设计
- 用户的问题是算法竞赛题或纯编程问题，不涉及深度学习建模
- 用户需要的是 MLOps 平台搭建或团队管理，而非具体技术方案

## 输入
```typescript
{
  task: {
    type: string                  // 任务类型(classification/detection/segmentation/generation/sequence-modeling/RL)
    domain: string                // 应用领域(CV/NLP/audio/multimodal/time-series)
    specificTask?: string         // 具体任务(如 image-classification → CIFAR-10)
  }
  data: {
    size?: string                 // 数据集规模(如 "10k 样本")
    quality?: string              // 数据质量描述(噪声、标注质量)
    distribution?: string         // 类别分布(平衡/长尾)
    labelAvailability?: string    // 标注情况(全监督/弱监督/无监督)
    augmentationNeeded?: boolean  // 是否需要数据增强
  }
  constraints: {
    computeResources?: string     // 计算资源(单GPU/多GPU/TPU/CPU)
    latencyRequirement?: string   // 延迟要求(如 "推理 <50ms")
    accuracyTarget?: string       // 精度目标(如 "top-1 acc >90%")
    modelSizeLimit?: string       // 模型大小限制(如 "<100MB")
    budget?: string               // 预算限制(影响预训练模型选择/计算资源)
  }
  challenges?: string[]           // 已知挑战(小样本/域迁移/实时性/长尾分布)
  existingSetup?: {
    framework?: string            // 现有框架(PyTorch/TensorFlow/JAX)
    baseModel?: string            // 已有基础模型
    trainingIssues?: string[]     // 当前训练问题
  }
}
```

## 输出
```typescript
{
  modelDesign: {
    architecture: string          // 推荐架构(CNN/Transformer/RNN/GNN/Diffusion/Hybrid)
    rationale: string             // 架构选择理由(归纳偏置/参数量/计算复杂度)
    pretrainedModel?: {
      name: string                // 预训练模型(如 ResNet-50, BERT-base, ViT-B/16)
      source: string              // 来源(HuggingFace/TorchVision/OpenAI)
      adaptation: string          // 适配方式(全量微调/LoRA/Adapter/冻结部分层)
    }
    tradeoff: string              // 从零训练 vs 微调 vs 迁移学习权衡分析
    codeFramework?: string        // PyTorch/TensorFlow 代码框架(如适用)
  }
  trainingStrategy: {
    dataProcessing: {
      normalization: string       // 归一化方法
      augmentation?: string[]     // 数据增强策略
      samplingStrategy?: string   // 采样策略(平衡采样/困难样本挖掘)
      batchSize: string           // 批大小建议(考虑内存和收敛)
    }
    optimizer: {
      type: string                // 优化器类型(SGD/AdamW/Lion)
      learningRate: string        // 初始学习率
      scheduler?: string          // 学习率调度(Cosine/StepLR/ReduceLROnPlateau)
      weightDecay?: string        // 权重衰减
    }
    trainingTechniques: {
      gradientClipping?: string   // 梯度裁剪阈值
      mixedPrecision?: boolean    // 混合精度训练(FP16/BF16)
      gradientAccumulation?: number // 梯度累积步数
      ema?: boolean               // 指数移动平均
    }
    regularization?: string[]     // 正则化方法(Dropout/Batch Norm/Label Smoothing)
  }
  debuggingGuidance: {
    visualization: string[]       // 可视化建议(loss曲线/权重分布/激活热图)
    checkpoints: string[]         // 检查点(权重范数/梯度范数/数据统计)
    commonIssues: {
      issue: string               // 问题类型(过拟合/欠拟合/NaN/不收敛)
      diagnosis: string           // 诊断方法
      solution: string[]          // 解决方案
    }[]
  }
  engineeringPractice: {
    modelCompression?: {
      pruning?: string            // 剪枝策略
      quantization?: string       // 量化方法(INT8/FP16)
      distillation?: string       // 知识蒸馏
    }
    inferenceOptimization?: {
      format: string[]            // 推理格式(ONNX/TensorRT/CoreML)
      batching?: string           // 批处理策略
      caching?: string            // 缓存策略
    }
    distributedTraining?: {
      strategy: string            // 分布式策略(DDP/DeepSpeed/FSDP)
      configuration: string       // 配置建议
    }
  }
  expectedPerformance?: {
    accuracy?: string             // 预期准确率
    trainingTime?: string         // 预期训练时间
    inferenceLatency?: string     // 预期推理延迟
  }
}
```

## 执行工作流

在开始执行前，复制以下清单，并在每一步完成后显式标记状态。

### Step 1: 任务需求分析
- 明确任务类型(分类/检测/分割/生成/序列建模/强化学习)
- 确认应用领域(CV/NLP/音频/多模态/时序)
- 量化数据情况(规模、质量、标注、分布)
- 识别约束条件(计算资源、延迟、精度、模型大小)
- 提取核心挑战(小样本、域迁移、长尾分布、实时性)

**反馈闭环**: 若关键信息缺失(如数据规模、计算资源)，必须询问用户补充，避免基于假设设计。

### Step 2: 模型架构选择
- 列出候选架构(CNN/Transformer/RNN/GNN/Diffusion/混合架构)
- 说明每个架构的归纳偏置(如 CNN 的局部性、Transformer 的长程依赖)
- 对比参数量、计算复杂度、数据需求
- 评估预训练模型可用性(HuggingFace/TorchVision/OpenAI)
- 权衡从零训练 vs 微调 vs 迁移学习

**架构选择原则**:
- **小数据(<10k)**: 优先迁移学习或微调预训练模型
- **中等数据(10k-100k)**: 微调 + 数据增强
- **大数据(>100k)**: 可考虑从零训练或轻度微调
- **实时性要求**: 优先轻量架构(MobileNet/EfficientNet/DistilBERT)

**反馈闭环**: 若用户对多个架构犹豫不决，提供 2-3 个方案的详细对比(精度/速度/资源消耗)。

### Step 3: 训练策略设计
- 数据处理: 归一化方法(ImageNet统计/自适应)、增强策略(AutoAugment/RandAugment)
- 优化器选择: SGD(视觉任务) vs AdamW(语言/多模态任务)
- 学习率调度: Warmup + Cosine Annealing(推荐) / StepLR / ReduceLROnPlateau
- 训练技巧: 梯度裁剪(防止爆炸)、混合精度(加速训练)、梯度累积(模拟大batch)
- 正则化: Dropout、权重衰减、标签平滑

**训练配方模板**:
```python
# 示例: 图像分类训练配方
optimizer = AdamW(model.parameters(), lr=1e-3, weight_decay=0.01)
scheduler = CosineAnnealingLR(optimizer, T_max=epochs)
scaler = GradScaler()  # 混合精度
# 数据增强: RandomCrop + RandomHorizontalFlip + ColorJitter
```

**反馈闭环**: 若用户报告训练不稳定(loss震荡/NaN)，立即检查学习率、数据归一化、梯度范数。

### Step 4: 实验执行与监控
- 建立基线实验(最简单可行的模型)
- 设置监控指标(训练/验证loss、准确率、梯度范数、学习率)
- 配置可视化工具(TensorBoard/WandB)
- 保存关键检查点(best model、last model、每N epoch)

**监控清单**:
- [ ] 训练loss平滑下降(无震荡)
- [ ] 验证loss与训练loss差距合理(<2x)
- [ ] 梯度范数稳定(无爆炸/消失)
- [ ] 学习率按预期衰减
- [ ] 数据加载无瓶颈(GPU利用率 >80%)

**反馈闭环**: 若验证loss早期不下降，检查数据pipeline(标签错误/归一化问题/增强过强)。

### Step 5: 问题诊断与调优
根据训练表现，诊断并解决常见问题:

#### 过拟合(训练acc高、验证acc低)
- **诊断**: 训练/验证loss曲线分离、验证准确率不再提升
- **解决方案**:
  - 增强数据增强(AutoAugment/MixUp/CutMix)
  - 增大 Dropout 率或权重衰减
  - 使用 Early Stopping
  - 减小模型容量(更少层/通道)
  - 获取更多训练数据

#### 欠拟合(训练/验证acc都低)
- **诊断**: 训练loss长期不下降、准确率低于随机猜测
- **解决方案**:
  - 增加模型容量(更深/更宽网络)
  - 降低正则化强度(减小 Dropout/权重衰减)
  - 检查数据质量(标注错误/类别不平衡)
  - 增加训练轮数或调整学习率

#### NaN/Inf(训练中断)
- **诊断**: loss突然变为NaN、梯度爆炸
- **解决方案**:
  - 降低学习率(减小10x)
  - 启用梯度裁剪(clip_grad_norm=1.0)
  - 检查输入数据(是否包含NaN/Inf)
  - 添加数值稳定性操作(如log时加epsilon)

#### 收敛慢
- **诊断**: loss缓慢下降、训练时间过长
- **解决方案**:
  - 增大学习率(谨慎,可能不稳定)
  - 使用更激进的学习率调度(如 OneCycleLR)
  - 启用混合精度训练(加速2-3x)
  - 优化数据加载(多进程、预取)

**反馈闭环**: 若多次调优仍无改善，重新审视 Step 2(架构选择)或 Step 1(任务定义)是否合理。

### Step 6: 模型压缩与部署优化(如适用)
- 评估压缩需求(模型大小/推理延迟是否满足约束)
- 选择压缩策略:
  - **剪枝**: 结构化剪枝(减少通道/层) vs 非结构化剪枝(减少参数)
  - **量化**: PTQ(训练后量化) vs QAT(量化感知训练), INT8/FP16
  - **蒸馏**: 使用大模型(教师)训练小模型(学生)
- 转换推理格式(ONNX/TensorRT/CoreML)
- 优化推理流程(批处理/算子融合/内存优化)

**压缩权衡**:
- **剪枝**: 精度下降中等(1-3%)、压缩比高(2-5x)
- **量化**: 精度下降小(<1%)、压缩比中等(2-4x)
- **蒸馏**: 精度下降可控、压缩比取决于学生模型大小

**反馈闭环**: 若压缩后精度下降超过预期(>5%)，考虑放宽压缩比或使用蒸馏补偿。

### Step 7: 性能验证与文档化
- 在测试集上验证最终性能(准确率、F1、mAP等)
- 测量推理延迟(单样本/批处理)
- 记录模型版本、训练配置、数据集版本
- 生成模型卡片(Model Card): 用途、性能、限制、伦理考虑

**性能基准对比**:
- 与公开基准(SOTA)对比
- 与简单基线(如线性模型)对比
- 分析失败案例(哪些类别/样本容易出错)

**反馈闭环**: 若性能未达目标且无明显改进空间，与用户讨论是否调整目标或增加资源投入。

## 失败处理

### 数据不足或质量差
- **现象**: 用户提供的数据量极小(<1000样本)或标注质量差
- **处理**: 
  - 推荐数据增强 + 预训练模型微调
  - 建议主动学习或半监督学习
  - 若预算允许,建议外包标注或合成数据生成

### 计算资源不足
- **现象**: 用户只有CPU或单个低端GPU,但任务需要大模型
- **处理**:
  - 推荐轻量架构(MobileNet/DistilBERT/TinyLlama)
  - 建议使用预训练模型 + 冻结部分层
  - 提供云平台方案(Colab/Kaggle/AWS SageMaker)

### 训练持续不收敛
- **现象**: 尝试多种学习率/优化器仍不收敛
- **处理**:
  - 检查数据pipeline(可视化batch样本)
  - 简化模型(从最小可行模型开始)
  - 检查标签是否正确(交叉验证标注质量)

### 部署后性能不达标
- **现象**: 模型在测试集表现好,但生产环境差
- **处理**:
  - 分析训练/生产数据分布差异(域偏移)
  - 使用域适应技术(DANN/CORAL)
  - 收集生产数据进行增量训练

### 过度依赖预训练模型
- **现象**: 用户无条件使用大预训练模型,忽略任务特性
- **处理**: 解释预训练模型适用场景,对比从零训练的成本/收益,推荐适中方案