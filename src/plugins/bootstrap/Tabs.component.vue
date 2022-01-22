<template>
  <ul class="nav nav-tabs">
    <li class="nav-item" v-for="name in tabs" :key="name">
      <strong><a
        class="nav-link"
        :class="current === name ? 'active' : ''"
        aria-current="page"
        href="#"
        @click.prevent="change(name)"
      >
        <slot :key="`${name}-title`"  :name="`${name}-title`"></slot>
      </a></strong>
    </li>
  </ul>
  <slot :key="current" :name="current"></slot>
</template>
<script>
export default {
  /* Usage:
    <Tabs>
    <template #tab1-title> Titulo da tab 1 </template>
    <template #tab1> conteudo da tab 1 </template>
    <template #tab2-title> Titulo da tab 2 </template>
    <template #tab2> conteudo da tab 2 </template>
  </Tabs>
  */
  data: () => ({ current: "", label: "" }),
  computed: {
    tabs: (t) => Object.keys(t.$slots).filter(t => !t.endsWith('-title')),
  },
  mounted() {
    this.current = this.tabs[0];
  },
  methods: {
    change(name) {
      this.current = name;
    },
  },
};
</script>