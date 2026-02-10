---
name: "Three.js着色器"
description: "Three.js 着色器：GLSL 基础、ShaderMaterial、uniforms、常用坐标空间与自定义效果。用户需要自定义视觉效果时调用。"
---

# Three.js 着色器

## 快速开始（ShaderMaterial）

```javascript
import * as THREE from "three";

const material = new THREE.ShaderMaterial({
  uniforms: {
    uTime: { value: 0 },
    uColor: { value: new THREE.Color(0x4aa3ff) },
  },
  vertexShader: `
    varying vec2 vUv;
    void main() {
      vUv = uv;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }
  `,
  fragmentShader: `
    uniform float uTime;
    uniform vec3 uColor;
    varying vec2 vUv;
    void main() {
      float wave = 0.5 + 0.5 * sin(uTime + vUv.x * 10.0);
      gl_FragColor = vec4(uColor * wave, 1.0);
    }
  `,
});

const mesh = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material);
scene.add(mesh);
```

## 核心概念

### 1) 常用矩阵与坐标空间
- `modelMatrix`：局部到世界
- `viewMatrix`：世界到相机
- `projectionMatrix`：相机到裁剪空间
- `modelViewMatrix`：`viewMatrix * modelMatrix`

### 2) uniforms / varyings
- `uniform`：CPU 传给 GPU 的参数（时间、颜色、贴图等）
- `varying`：顶点着色器传给片元着色器的插值数据（如 UV、法线等）

### 3) 与 Three.js 管线协作
- 常见需求：雾、阴影、色彩空间、灯光
- 如果需要完整灯光/PBR，一般考虑 `onBeforeCompile` 或使用 NodeMaterial/后处理管线

## 常用模式

### 1) 时间驱动动画

```javascript
const clock = new THREE.Clock();
function animate() {
  requestAnimationFrame(animate);
  material.uniforms.uTime.value += clock.getDelta();
  renderer.render(scene, camera);
}
```

### 2) 纹理采样
- 在 fragment shader 中用 `sampler2D` + `texture2D/texture`

## 性能提示
- shader 分支与循环会增加成本，尽量保持简单
- 纹理采样次数是热点，能少就少
- 后处理效果往往更容易堆成本，注意分辨率与 pass 数量

## 相关技能
- Three.js 后处理
- Three.js 材质
- Three.js 纹理
