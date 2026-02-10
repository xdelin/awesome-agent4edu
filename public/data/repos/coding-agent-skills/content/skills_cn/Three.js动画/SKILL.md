---
name: "Three.js动画"
description: "Three.js 动画：AnimationMixer、关键帧、骨骼/蒙皮、morph targets 与混合。用户需要播放或混合模型动画时调用。"
---

# Three.js 动画

## 快速开始（AnimationMixer）

```javascript
import * as THREE from "three";

const clock = new THREE.Clock();

const mixer = new THREE.AnimationMixer(gltf.scene);
const clip = gltf.animations[0];
const action = mixer.clipAction(clip);
action.play();

function animate() {
  requestAnimationFrame(animate);
  mixer.update(clock.getDelta());
  renderer.render(scene, camera);
}
animate();
```

## 核心概念

### 1) Clip / Action / Mixer
- `AnimationClip`：动画数据（轨道集合）
- `AnimationAction`：某个 clip 在 mixer 上的播放实例（可暂停/淡入淡出/权重）
- `AnimationMixer`：驱动更新与混合

### 2) 循环与时间缩放

```javascript
action.setLoop(THREE.LoopRepeat, Infinity);
action.timeScale = 1.0;
```

### 3) 动画混合（淡入淡出）

```javascript
actionA.reset().fadeIn(0.2).play();
actionB.fadeOut(0.2);
```

## 常用模式

### 1) 动画状态机
- 用状态机管理 “idle/run/jump”等动作，统一处理切换与淡入淡出

### 2) 多个动画并行
- 同一个 mixer 可以同时运行多个 action，通过权重实现叠加/混合

## 性能提示
- `mixer.update(dt)` 建议使用稳定 `dt`（Clock.getDelta）
- 对于大量角色，注意骨骼数量与动画采样成本
- 不需要的 action 及时 `stop()`，并避免每帧创建 action

## 相关技能
- Three.js 加载器
- Three.js 基础
