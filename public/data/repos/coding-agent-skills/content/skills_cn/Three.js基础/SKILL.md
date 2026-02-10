---
name: "Three.js基础"
description: "Three.js 基础能力：场景/相机/渲染器、坐标系、Object3D 层级与变换。用户需要搭建 3D 场景或理解基础概念时调用。"
---

# Three.js 基础

## 快速开始

```javascript
import * as THREE from "three";

const scene = new THREE.Scene();

const camera = new THREE.PerspectiveCamera(
  75,
  window.innerWidth / window.innerHeight,
  0.1,
  1000,
);
camera.position.set(0, 0, 5);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
document.body.appendChild(renderer.domElement);

const cube = new THREE.Mesh(
  new THREE.BoxGeometry(1, 1, 1),
  new THREE.MeshStandardMaterial({ color: 0x00ff00 }),
);
scene.add(cube);

scene.add(new THREE.AmbientLight(0xffffff, 0.5));
const dirLight = new THREE.DirectionalLight(0xffffff, 1);
dirLight.position.set(5, 5, 5);
scene.add(dirLight);

function animate() {
  requestAnimationFrame(animate);
  cube.rotation.y += 0.01;
  renderer.render(scene, camera);
}
animate();

window.addEventListener("resize", () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});
```

## 核心概念

### 1) 场景、相机、渲染器
- `Scene`：容器，持有所有可渲染对象/灯光/雾等
- `Camera`：决定视角与投影（透视/正交）
- `WebGLRenderer`：将场景+相机渲染到 `<canvas>`

### 2) Object3D 层级
- 大部分对象都继承自 `Object3D`（如 `Mesh`、`Group`、`Light`）
- 父子层级会影响变换：子对象的世界矩阵由父矩阵叠加

```javascript
const group = new THREE.Group();
group.position.set(0, 1, 0);
scene.add(group);

const child = new THREE.Mesh(new THREE.BoxGeometry(), new THREE.MeshBasicMaterial());
child.position.set(2, 0, 0);
group.add(child);
```

### 3) 变换与坐标系
- 局部变换：`position` / `rotation` / `scale`
- 需要时用 `object.updateMatrixWorld()` 刷新世界矩阵
- 常用辅助：`AxesHelper`、`GridHelper`

## 常用模式

### 1) 渲染循环与时间步进
- 用 `THREE.Clock` 或 `performance.now()` 计算 `deltaTime`，让动画与帧率解耦

```javascript
const clock = new THREE.Clock();
function animate() {
  requestAnimationFrame(animate);
  const dt = clock.getDelta();
  cube.rotation.y += dt;
  renderer.render(scene, camera);
}
```

### 2) 资源释放（避免显存泄漏）
- 几何体/材质/纹理等需要在不再使用时 `dispose()`

```javascript
geometry.dispose();
material.dispose();
texture.dispose();
```

## 性能提示
- 像素比：`renderer.setPixelRatio(Math.min(devicePixelRatio, 2))` 通常更稳
- 避免每帧创建新对象（如 `new Vector3()`），尽量复用
- 大量对象渲染优先考虑 Instancing（见“Three.js 几何体”）

## 相关技能
- Three.js 几何体
- Three.js 材质
- Three.js 灯光
- Three.js 纹理
- Three.js 交互
