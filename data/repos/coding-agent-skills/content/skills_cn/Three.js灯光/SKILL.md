---
name: "Three.js灯光"
description: "Three.js 灯光：灯光类型、阴影、环境光照与常用调试辅助。用户需要布光或处理阴影问题时调用。"
---

# Three.js 灯光

## 快速开始

```javascript
import * as THREE from "three";

scene.add(new THREE.AmbientLight(0xffffff, 0.3));

const dir = new THREE.DirectionalLight(0xffffff, 1);
dir.position.set(5, 8, 5);
dir.castShadow = true;
scene.add(dir);
```

## 核心概念

### 1) 常用灯光类型
- `AmbientLight`：全局环境光（无方向）
- `DirectionalLight`：平行光（类似太阳），最常用于阴影
- `PointLight`：点光源（近似灯泡）
- `SpotLight`：聚光灯（锥形）
- `HemisphereLight`：天空/地面半球光（适合户外氛围）

### 2) 阴影基础
- 需要同时满足：
  - `renderer.shadowMap.enabled = true`
  - 光源 `castShadow = true`
  - 投射物体 `castShadow = true`
  - 接收物体 `receiveShadow = true`

```javascript
renderer.shadowMap.enabled = true;
mesh.castShadow = true;
ground.receiveShadow = true;
```

### 3) 阴影质量与范围
- 阴影模糊/锯齿通常来自：
  - 阴影贴图分辨率 `light.shadow.mapSize`
  - 阴影相机范围（DirectionalLight 的 orthographic shadow camera）

## 常用模式

### 1) 使用 helper 调试

```javascript
scene.add(new THREE.DirectionalLightHelper(dir, 0.5));
scene.add(new THREE.CameraHelper(dir.shadow.camera));
```

### 2) 环境光照与 HDRI
- PBR 场景常用环境贴图作为 IBL（见“Three.js 纹理”）

## 性能提示
- 阴影开销大：减少投射阴影的灯光数量与投射物体数量
- 缩小阴影相机范围可显著提升质量与稳定性
- 优先使用 DirectionalLight + 单张阴影贴图满足大多数场景

## 相关技能
- Three.js 材质
- Three.js 纹理
- Three.js 基础
