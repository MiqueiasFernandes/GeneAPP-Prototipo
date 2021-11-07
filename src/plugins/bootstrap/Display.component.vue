<template>
  <div v-if="lead">
    <p class="lead">
      <Icon :name="ico" class="icon-lead" v-if="ico" /><slot></slot>
    </p>
  </div>
  <div v-else>
    <h1 :class="classObj">
      <Icon
        :class="subtitle ? 'icon-subtitle' : 'icon-title'"
        :name="ico"
        v-if="ico"
      /><slot></slot>
    </h1>
  </div>
</template>
<script>
export default {
  props: {
    lead: Boolean,
    subtitle: Boolean,
    ico: String,
    size: Number
  },
  computed: {
    classObj: (t) =>
      `display-${
        t.size ? t.size : t.subtitle ? t.$bootstrap_conf.subtitle : t.$bootstrap_conf.title
      }`,
  },
};
</script>
<style scoped>
.icon-title {
  position: absolute;
  transform: rotate(-25deg) scale(4);
  opacity: 0.08;
  left: 10rem;
  top: 8rem;
}
.icon-subtitle {
  width: 3rem;
  height: 3rem;
  padding-bottom: 0.6rem;
}
.icon-lead {
  margin: 0 0.5rem 0.1rem 0;
  width: 1.5rem;
  height: 1.5rem;
}
</style>