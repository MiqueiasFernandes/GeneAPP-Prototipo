<template>
  <svg
    :innerHTML="icon"
    :class="color ? ['bi', 'text-' + color] : 'bi'"
    :width="sm ? '20' : lg ? '38' : size"
    :height="sm ? '20' : lg ? '38' : size"
    fill="currentColor"
  ></svg>
</template>

<script>
export default {
  props: {
    fill: Boolean,
    color: {
      type: String,
    },
    name: {
      type: String,
    },
    size: {
      type: String,
      default: "32",
    },
    sm: Boolean,
    lg: Boolean,
  },
  data: () => ({ icon: "" }),
  beforeMount() {
    this.setIcon();
  },
  watch: {
    name() {
      this.setIcon();
    },
  },
  methods: {
    setIcon() {
      this.$bootstrap_icons.then(
        (icons) =>
          (this.icon = icons[`${this.name}${this.fill ? "-fill" : ""}`])
      );
    },
  },
};
</script>