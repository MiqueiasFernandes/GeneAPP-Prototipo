<template>
  <span :class="classObj">
    <slot></slot>
  </span>
</template>
<script>
export default {
  props: {
    grow: Boolean,
    sm: Boolean,
    lg: Boolean,
    color: String,
  },
  computed: {
    classObj: (t) => {
      const type = `spinner-${t.grow ? "grow" : "border"}`;
      const classes = [type];
      if (t.sm) {
        classes.push(`${type}-sm`);
      } else if (t.lg) {
        classes.push("spiner-lg");
      }
      if (t.color) {
        classes.push(`text-${t.color}`);
      }
      return classes;
    },
  },
};
</script>
<style scoped>
.spiner-lg {
  width: 3rem;
  height: 3rem;
}
</style>