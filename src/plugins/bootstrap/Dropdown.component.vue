<template>
  <div class="dropdown">
    <button
      class="btn dropdown-toggle"
      :class="`btn-${color}`"
      type="button"
      id="dropdownMenuButton"
      data-toggle="dropdown"
      aria-expanded="false"
    >
      <slot>{{ title }}</slot>
    </button>
    <ul class="dropdown-menu" aria-labelledby="dropdownMenuButton">
      <li v-for="(option, i) in options" v-bind:key="i">
        <a class="dropdown-item" @click="send(option)">
          {{ option.label || option }}
        </a>
      </li>
    </ul>
  </div>
</template>

<script>
export default {
  emits: ["choose"],
  props: {
    title: String,
    options: Array,
    color: { type: String, default: "primary" },
  },
  methods: {
    send(option) {
      if (option.click) {
        option.click();
      } else {
        this.$emit("choose", option.id || option);
      }
    },
  },
};
</script>