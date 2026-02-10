---
name: "Three.js几何体"
description: "Three.js 几何体：内置几何、BufferGeometry、自定义属性、Instancing。用户需要创建/优化网格几何时调用。"
---

# Three.js 几何体

## 快速开始

```javascript
import * as THREE from "three";

const geometry = new THREE.BoxGeometry(1, 1, 1);
const material = new THREE.MeshStandardMaterial({ color: 0x4aa3ff });
const mesh = new THREE.Mesh(geometry, material);
scene.add(mesh);
```

## 核心概念

### 1) Geometry vs BufferGeometry
- 现代 Three.js 主要使用 `BufferGeometry`
- 顶点数据存储在 `geometry.attributes` 中（如 `position`、`normal`、`uv`）

### 2) 常见属性
- `position`：每个顶点的 xyz（`Float32BufferAttribute`）
- `normal`：法线（影响光照）
- `uv`：纹理坐标
- `index`：索引（复用顶点）

```javascript
const positions = new Float32Array([
  0, 0, 0,
  1, 0, 0,
  0, 1, 0,
]);
const geom = new THREE.BufferGeometry();
geom.setAttribute("position", new THREE.BufferAttribute(positions, 3));
geom.computeVertexNormals();
```

### 3) Instancing（大批量同材质对象）
- 用 `InstancedMesh` 将大量相同几何+材质的对象合批渲染

```javascript
const count = 1000;
const g = new THREE.SphereGeometry(0.05, 16, 16);
const m = new THREE.MeshStandardMaterial({ color: 0xffffff });
const instanced = new THREE.InstancedMesh(g, m, count);
const mat4 = new THREE.Matrix4();
for (let i = 0; i < count; i++) {
  mat4.makeTranslation(Math.random() * 10 - 5, Math.random() * 2, Math.random() * 10 - 5);
  instanced.setMatrixAt(i, mat4);
}
instanced.instanceMatrix.needsUpdate = true;
scene.add(instanced);
```

## 常用模式

### 1) 更新几何属性
- 修改 attribute 后需要设置 `needsUpdate = true`

```javascript
const pos = geom.attributes.position;
pos.setXYZ(0, 0, 0, 0.2);
pos.needsUpdate = true;
```

### 2) 计算包围盒/包围球
- 用于视锥裁剪、碰撞、相机 framing 等

```javascript
geom.computeBoundingBox();
geom.computeBoundingSphere();
```

## 性能提示
- 能用 `index` 就用索引：减少顶点重复
- 大量重复对象优先 Instancing
- 动态几何尽量减少顶点数量与更新频率

## 相关技能
- Three.js 基础
- Three.js 材质
- Three.js 灯光
- Three.js 纹理
