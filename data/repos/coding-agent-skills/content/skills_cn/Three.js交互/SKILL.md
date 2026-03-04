---
name: "Three.js交互"
description: "Three.js 交互：Raycaster、鼠标/触摸输入、对象拾取、高亮与控制器集成。用户需要点击/拖拽/选择等交互时调用。"
---

# Three.js 交互

## 快速开始（Raycaster 拾取）

```javascript
import * as THREE from "three";

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();

function onPointerMove(e) {
  pointer.x = (e.clientX / window.innerWidth) * 2 - 1;
  pointer.y = -(e.clientY / window.innerHeight) * 2 + 1;
}
window.addEventListener("pointermove", onPointerMove);

function animate() {
  requestAnimationFrame(animate);

  raycaster.setFromCamera(pointer, camera);
  const hits = raycaster.intersectObjects(scene.children, true);
  if (hits.length > 0) {
    const hit = hits[0].object;
    // 在这里做高亮/提示等交互反馈
  }

  renderer.render(scene, camera);
}
animate();
```

## 核心概念

### 1) Raycaster
- `setFromCamera(ndc, camera)`：将屏幕坐标（NDC）投射到 3D 射线
- `intersectObjects(objs, recursive)`：返回相交信息（距离、uv、face 等）

### 2) 事件与坐标转换
- 鼠标/触摸/指针事件 → 屏幕像素坐标 → NDC（[-1, 1]）

### 3) 常见交互反馈
- 悬停高亮（改变材质/outline pass）
- 点击选中（保存选中对象引用）
- 拖拽与 gizmo（配合控制器或自定义平面约束）

## 常用模式

### 1) 仅对可交互对象做检测
- 建议维护一个 `pickables` 数组，避免每次对整个场景树检测

```javascript
const pickables = [meshA, meshB];
const hits = raycaster.intersectObjects(pickables, true);
```

### 2) 与控制器共存
- 若使用 OrbitControls/PointerLockControls 等，注意事件优先级与冲突处理

## 性能提示
- Raycaster 每帧检测成本随对象数量增长；尽量减少检测集合
- 对于大量实例对象，考虑使用更轻量的拾取策略（如 GPU picking 或 BVH）

## 相关技能
- Three.js 基础
- Three.js 几何体
- Three.js 后处理
