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
        {{ name }}
      </a></strong>
    </li>
  </ul>
  <slot :key="current" :name="current"></slot>
</template>
<script>
export default {
  data: () => ({ current: "" }),
  computed: {
    tabs: (t) => Object.keys(t.$slots),
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