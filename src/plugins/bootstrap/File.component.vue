<template>
  <div class="form-file" :class="{ 'form-file-lg': lg }">
    <label class="form-file-label" for="fileInput"
      v-if="opts">
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
    <input
      type="file"
      class="form-control"
      id="fileInput"
      :disabled="disabled"
      :multiple="max > 1"
      ref="file"
      v-on:change="setFile"
    />
    <Button
      v-if="opts"
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
      v-if="opts && clear"
      >clear</Button
    >
  </div>
</template>
<script>
export default {
  emits: ["files", "parsed", "load"],
  props: {
    opts: {
      type: Boolean,
      default: true,
    },
    disabled: Boolean,
    lg: Boolean,
    img: Boolean,
    ico: String,
    txt: Boolean,
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
        if (this.txt) {
          this.toText();
        }
        this.removed = !this.files;
      }
    },

    validFile(or = null) {
      const files = [...this.$refs.file.files];
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
          this.extension &&
          this.extension.length > 0 &&
          this.extension.every((t) => !file.name.endsWith(t))
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

    toText() {
      const parsed = this.files.map((f) => ({
        file: f,
        status: 0,
        txt: null,
      }));
      const int = setInterval(() => {
        if (parsed.every((f) => f.status > 1)) {
          clearInterval(int);
          this.$emit("parsed", parsed);
        } else {
          if (!parsed.some((f) => f.status === 1)) {
            const file = parsed.filter((f) => f.status < 1)[0];
            file.status++;
            var reader = new FileReader();
            reader.onload = (reader) => {
              file.txt = reader.target.result;
              file.status++;
              this.$emit("load", [
                parsed.filter((p) => p.status > 1).length,
                parsed.length,
              ]);
            };
            reader.readAsText(file.file);
          }
        }
      }, 300);
    },
  },
};
</script>