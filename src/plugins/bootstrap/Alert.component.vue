<template>
  <div
    class="alert alert-dismissible fade"
    :class="alert_color"
    role="alert"
    ref="alert"
    :hidden="!display"
  >
    <slot></slot>
    <button
    v-if="!noclose"
      type="button"
      class="btn-close"
      aria-label="Close"
      @click="hide"
    ></button>
  </div>
</template>
<script>
export default {
  emits: ['close'],
  props: {
    color: {
      type: String,
      default: "warning",
    },
    success: Boolean,
    primary: Boolean,
    danger: Boolean,
    hidden: Boolean,
    noclose: Boolean
  },
  computed: {
    alert_color: (t) => {
      const color = t.success
        ? "success"
        : t.primary
        ? "primary"
        : t.danger
        ? "danger"
        : t.color;
      return [`alert-${color}`, t.display ? "show" : ""];
    },
    instance: (t) => t.$bootstrap.Alert.getInstance(t.$refs.alert),
  },
  data: () => ({ display: false }),
  mounted() {
    new this.$bootstrap.Alert(this.$refs.alert);
    if (!this.hidden) {
      this.show();
    }
  },
  methods: {
    show() {
      this.display = true;
    },
    hide() {
      if (this.hidden) {
        this.display = false;
      } else {
        this.instance.close();
      }
      this.$emit('close')
    },
  },
};
</script>