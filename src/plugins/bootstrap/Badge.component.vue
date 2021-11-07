<template>
  <span class="badge" :class="classObj">
    <slot></slot>
    <button
      :disabled="disabled"
      @click="close"
      v-if="closeable"
      type="button"
      class="btn-close"
      :class="is_dark ? 'btn-close-white' : ''"
      aria-label="Close"
    ></button>
  </span>
</template>
<script>
export default {
  emits: ["close"],
  props: {
    color: {
      type: String,
      default: "secondary",
    },
    danger: Boolean,
    warning: Boolean,
    success: Boolean,
    round: Boolean,
    closeable: Boolean,
    disabled: Boolean,
  },
  computed: {
    parsed_color: (t) =>
      `bg-${
        t.danger
          ? "danger"
          : t.warning
          ? "warning"
          : t.success
          ? "success"
          : t.color
      }`,
    classObj: (t) => [
      t.parsed_color,
      t.round ? "rounded-pill" : "",
      t.disabled ? "disabled" : "",
    ],
    is_dark: (t) =>
      ["secondary", "dark"].includes(t.parsed_color.split("-")[1]),
  },
  methods: {
    close() {
      this.$emit("close");
    },
  },
};
</script>
<style scoped>
.disabled {
  opacity: .8;
}
</style>