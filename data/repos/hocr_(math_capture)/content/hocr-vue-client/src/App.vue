<template>
  <div class="container">
    <div class="toolbar">
      <label class="upload-btn">
        Upload
        <input type="file" accept="image/png, image/jpeg" @change="onFileChange" />
      </label>
      <input v-model="instruction" class="prompt-input" placeholder="add extra prompt..." />
      <button class="check-btn" @click="onSubmit">check</button>
    </div>

    <div class="grid">
      <div class="preview">
        <div v-if="previewUrl">
          <img :src="previewUrl" alt="uploaded image" />
        </div>
        <div v-else class="placeholder">uploaded image</div>
      </div>

      <div class="results">
        <div class="extracted">
          <h3>VLM answer result</h3>
          <pre>{{ result }}</pre>
        </div>
        <div class="extracted">
          <h3>Extracted latex code</h3>
          <pre>{{ latex }}</pre>
        </div>
        <div class="rendered">
          <h3>Rendered LaTeX</h3>
          <div v-html="renderedHtml" class="math"></div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import axios from 'axios';
import katex from 'katex';

export default {
  data() {
    return {
      file: null,
      previewUrl: '',
      instruction: '',
      latex: ''
    };
  },
  computed: {
    renderedHtml() {
      try {
        return this.latex
          ? katex.renderToString(this.latex, { throwOnError: false })
          : '';
      } catch (e) {
        return '<span style="color: red">渲染错误</span>';
      }
    }
  },
  methods: {
    onFileChange(e) {
      const f = e.target.files[0];
      if (!f) return;
      this.file = f;
      this.previewUrl = URL.createObjectURL(f);
    },
    async onSubmit() {
      if (!this.file) {
        alert('请选择图片');
        return;
      }
      const form = new FormData();
      form.append('file', this.file);
      form.append('instruction', this.instruction);

      try {
        const resp = await axios.post('/upload', form);
        this.result = resp.data.result || '';
        this.latex = resp.data.latex || '';
      } catch (err) {
        console.error(err);
        alert('识别失败');
      }
    }
  }
};
</script>

<style scoped>
.container { padding: 20px; font-family: Arial, sans-serif; }
.toolbar { display: flex; align-items: center; gap: 10px; background: #eee; padding: 10px; }
.upload-btn { position: relative; background: #fff; border: 1px solid #ccc; padding: 6px 12px; cursor: pointer; }
.upload-btn input[type=file] { position: absolute; left: 0; top: 0; width: 100%; height: 100%; opacity: 0; cursor: pointer; }
.prompt-input { flex: 1; padding: 6px; border: 1px solid #ccc; }
.check-btn { padding: 6px 12px; border: 1px solid #007bff; background: #007bff; color: #fff; cursor: pointer; }
.grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 20px; }
.preview, .extracted, .rendered { background: #f5f5f5; padding: 20px; border: 1px solid #ddd; min-height: 200px; }
.placeholder { color: #999; text-align: center; line-height: 160px; }
.extracted pre { white-space: pre-wrap; word-break: break-word; }
.math { margin-top: 10px; }
</style>