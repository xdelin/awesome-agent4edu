---
name: "Three.js加载器"
description: "Three.js 加载器：GLTF/GLB、纹理加载、异步与缓存、常见坑（色彩空间/坐标系/资源路径）。用户需要加载模型或资源时调用。"
---

# Three.js 加载器

## 快速开始（GLTFLoader）

```javascript
import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js";

const loader = new GLTFLoader();
loader.load(
  "/model.glb",
  (gltf) => {
    scene.add(gltf.scene);
  },
  undefined,
  (err) => {
    console.error(err);
  },
);
```

## 核心概念

### 1) 常用加载器
- `GLTFLoader`：加载 glTF/glb（最常用的 3D 资产格式）
- `TextureLoader`：加载图片纹理

### 2) Draco 压缩（可选）
- 需要 `DRACOLoader` 配合，并配置解码器路径

### 3) 加载后的处理
- 模型通常需要：
  - 设置缩放/位置/旋转
  - 遍历场景树处理材质与阴影

```javascript
gltf.scene.traverse((obj) => {
  if (obj.isMesh) {
    obj.castShadow = true;
    obj.receiveShadow = true;
  }
});
```

## 常用模式

### 1) Promise 封装
- 统一加载流程与错误处理，便于并行加载与进度管理

### 2) 资源缓存
- 相同 URL 的模型/纹理尽量复用，避免重复下载与重复解析

## 性能提示
- 模型面数与材质数量决定渲染开销，必要时做简化与合批
- 优先使用 glb（二进制）与压缩纹理（KTX2）
- 加载完成后再启用后处理/高质量阴影，缩短首屏时间

## 相关技能
- Three.js 纹理
- Three.js 动画
- Three.js 后处理
