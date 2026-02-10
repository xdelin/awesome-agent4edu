---
name: "Three.js材质"
description: "Three.js 材质：PBR（Standard/Physical）、基础材质、透明度、贴图与渲染状态。用户需要选择/调试材质效果时调用。"
---

# Three.js 材质

## 快速开始

```javascript
import * as THREE from "three";

const material = new THREE.MeshStandardMaterial({
  color: 0xffffff,
  metalness: 0.2,
  roughness: 0.6,
});
const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 32), material);
scene.add(mesh);
```

## 核心概念

### 1) 常用材质选型
- `MeshStandardMaterial`：通用 PBR（推荐默认）
- `MeshPhysicalMaterial`：更高级 PBR（clearcoat、transmission 等）
- `MeshBasicMaterial`：不受光照影响（调试/卡通 UI 物体）
- `MeshPhongMaterial`：传统高光模型（非 PBR）

### 2) 透明与混合
- 透明一般需要：
  - `material.transparent = true`
  - `material.opacity = 0 ~ 1`
- 透明物体可能出现排序问题，可按需调整：
  - `material.depthWrite = false`
  - `renderOrder`

### 3) 双面与背面剔除
- 默认开启背面剔除（性能更好）
- 需要双面渲染时：

```javascript
material.side = THREE.DoubleSide;
```

### 4) 贴图与色彩空间
- 颜色贴图通常需要 sRGB 色彩空间

```javascript
const loader = new THREE.TextureLoader();
const map = loader.load("/albedo.jpg");
map.colorSpace = THREE.SRGBColorSpace;
material.map = map;
material.needsUpdate = true;
```

## 常用模式

### 1) 环境光照（IBL）
- PBR 材质通常需要 `scene.environment` 提供环境反射/漫反射（HDRI）

```javascript
scene.environment = envMap;
```

### 2) 统一调参
- 将材质参数集中管理，方便 GUI/调试（如 roughness/metalness）

## 性能提示
- 材质种类越多，越难合批；尽量复用材质实例
- 透明物体与后处理成本高，慎用大面积透明
- 贴图分辨率与数量会直接影响显存与带宽

## 相关技能
- Three.js 纹理
- Three.js 灯光
- Three.js 着色器
