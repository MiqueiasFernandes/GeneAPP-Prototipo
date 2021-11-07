<template>
  <div class="form-file" :class="{ 'form-file-lg': lg }">
    <input
      type="file"
      class="form-file-input"
      id="fileInput"
      :disabled="disabled"
      :multiple="max > 1"
      ref="file"
      v-on:change="setFile"
    />
    <label class="form-file-label" for="fileInput">
      <span class="form-file-text">
        <Badge v-if="files && files.length > 1" secondary round>
          {{ files.length }}
        </Badge>
        {{ text }}</span
      >
      <span class="form-file-button">
        <Icon :name="ico ? ico : img ? 'camera' : 'paperclip'" sm />
      </span>
    </label>
    <Button
      type="button"
      sm
      color="none"
      class="mt-1"
      @click="restore"
      ico="arrow-clockwise"
      :disabled="!fileWasChanged"
      >reset</Button
    >
    <Button
      type="button"
      sm
      color="none"
      class="mt-1"
      @click="remove"
      ico="x"
      v-if="clear"
      >clear</Button
    >
  </div>
</template>
<script>
export default {
  emits: ["files"],
  props: {
    disabled: Boolean,
    lg: Boolean,
    img: Boolean,
    ico: String,
    min: {
      type: Number,
      default: 0,
    },
    max: {
      type: Number,
      default: 1,
    },
    size: {
      type: Number,
      default: 1000,
    },
    types: Array,
    clear: {
      type: Boolean,
      default: true,
    },
    value: {
      type: String,
      default: "Choose file...",
    },
  },
  computed: {
    extension: (t) => (t.types ? t.types : t.img ? ["image/"] : null),
    text: (t) =>
      t.removed
        ? "-"
        : t.files && t.files.length > 0
        ? Array(...t.files)
            .map((f) => f.name)
            .join(", ")
        : t.value,
    fileWasRemoved: (t) => t.removed,
    fileWasChanged: (t) => (t.files && t.files.length > 0) || t.removed,
    selectedFile: (t) => t.validFile([null])[0],
    selectedFiles: (t) => t.validFile([]),
  },
  data: () => ({ files: [], removed: false }),
  methods: {
    restore() {
      this.removed = false;
      this.resetInput();
      this.files = null;
    },

    remove() {
      this.files = null;
      this.resetInput();
      this.removed = true;
    },

    setFile(event) {
      if (event.target && event.target.files) {
        this.files = this.validFile();
        this.removed = !this.files;
      }
    },

    validFile(or = null) {
      const files = this.$refs.file.files;
      if (!files) {
        return or;
      }
      if (files.length < this.min) {
        this.resetInput();
        alert("selecione no minimo " + this.min + " arquivo");
        return or;
      }
      if (files.length > this.max) {
        this.resetInput();
        alert("selecione no maximo " + this.max + " arquivo");
        return or;
      }
      let reset = false;
      const erros = [];

      files.forEach((file) => {
        if (
          !file.type ||
          (this.extension &&
            this.extension.every((t) => file.type.indexOf(`${t}`) !== 0))
        ) {
          if (!erros.includes(file.type)) {
            erros.push(file.type);
            alert("Tipo invalido de arquivo " + file.type);
          }
          reset = true;
          return;
        }
        if (file.size > 1000 * this.size) {
          if (!erros.includes(this.size)) {
            const size = this.size < 1024 ? this.size : this.size / 1024;
            const suffix = this.size < 1024 ? "Kb" : "Mb";
            alert("Selecione arquivo menor que " + size + suffix);
            erros.push(this.size);
          }
          reset = true;
        }
      });

      if (reset) {
        this.resetInput();
      } else {
        this.$emit("files", files);
        return files;
      }
      return or;
    },

    resetInput() {
      this.$refs.file.files = null;
      this.$refs.file.value = null;
    },
  },
};
</script>