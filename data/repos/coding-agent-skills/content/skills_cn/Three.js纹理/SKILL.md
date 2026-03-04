---
name: "Three.js纹理"
description: "Three.js 纹理：TextureLoader、UV、wrap/repeat、色彩空间、环境贴图与渲染目标。用户需要加载/调试贴图时调用。"
---

# Three.js 纹理

## 快速开始

```javascript
import * as THREE from "three";

const loader = new THREE.TextureLoader();
const map = loader.load("/albedo.jpg");
map.colorSpace = THREE.SRGBColorSpace;

const material = new THREE.MeshStandardMaterial({ map });
```

## 核心概念

### 1) wrap / repeat / offset

```javascript
map.wrapS = THREE.RepeatWrapping;
map.wrapT = THREE.RepeatWrapping;
map.repeat.set(2, 2);
map.offset.set(0.1, 0.2);
map.needsUpdate = true;
```

### 2) 过滤与 mipmaps
- `minFilter` / `magFilter` 影响缩放时清晰度与性能
- 贴图分辨率越高，显存占用越大

### 3) 色彩空间与编码
- 颜色贴图（albedo/baseColor）通常用 `SRGBColorSpace`
- 非颜色数据（normal/roughness/metalness/ao）一般保持线性（不要设 sRGB）

### 4) 环境贴图（envMap）
- PBR 材质非常依赖环境贴图（反射与间接光）
- 通常将环境贴图赋给 `scene.environment` 与 `scene.background`

## 常用模式

### 1) 贴图缓存与复用
- 同一 URL 的贴图尽量只加载一次，避免重复占用显存

### 2) RenderTarget（渲染到纹理）
- 用于后处理、反射、离屏渲染等高级效果（见“Three.js 后处理”）

## 性能提示
- 控制贴图尺寸与数量，尽量使用压缩纹理（KTX2 等）
- 小纹理/像素风可使用 `NearestFilter`，但注意闪烁与锯齿
- 大场景建议用 HDRI + PMREM（环境预滤）提升效果与一致性

## 相关技能
- Three.js 材质
- Three.js 灯光
- Three.js 后处理
