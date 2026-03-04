---
name: "Three.js后处理"
description: "Three.js 后处理：EffectComposer、常见 Pass（Bloom/DOF 等）、屏幕特效与自定义 Pass。用户需要屏幕后处理链路时调用。"
---

# Three.js 后处理

## 快速开始（EffectComposer）

```javascript
import * as THREE from "three";
import { EffectComposer } from "three/addons/postprocessing/EffectComposer.js";
import { RenderPass } from "three/addons/postprocessing/RenderPass.js";

const composer = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, camera));

function animate() {
  requestAnimationFrame(animate);
  composer.render();
}
animate();

window.addEventListener("resize", () => {
  composer.setSize(window.innerWidth, window.innerHeight);
});
```

## 核心概念

### 1) 为什么用后处理
- 将渲染结果作为纹理再处理，可实现 Bloom、色调映射、景深、抗锯齿等屏幕特效

### 2) Pass 链
- `RenderPass` 通常是第一步
- 后续叠加各种 pass（注意顺序会影响效果）

## 常用模式

### 1) 控制后处理分辨率
- 在性能敏感场景可降低 composer 分辨率（特别是移动端）

### 2) 与透明/深度配合
- 部分效果需要 depth texture 或正确的 alpha/深度写入设置

## 性能提示
- Pass 越多成本越高；优先合并效果或减少高分辨率 pass
- Bloom/DOF 等往往涉及多次模糊与采样，注意参数与分辨率

## 相关技能
- Three.js 着色器
- Three.js 纹理
