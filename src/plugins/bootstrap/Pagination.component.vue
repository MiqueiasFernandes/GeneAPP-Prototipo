<template>
  <nav aria-label="pages-navigation">
    <ul class="pagination pagination-sm justify-content-end" v-if="pages > 1">
      <li
        v-for="(page, i) in parsed_pages"
        v-bind:key="i"
        class="page-item"
        :class="{ active: page.current }"
        :aria-current="page.current ? 'page' : ''"
      >
        <a class="page-link active" v-if="page.current" href="#">
          {{ page.label }}
          <!-- <span class="visually-hidden">(current)</span> -->
        </a>
        <a v-else  class="page-link" @click.prevent="$emit('page', page.label)"  href="#">{{
            page.label
          }}</a>
      </li>
    </ul>
  </nav>
</template>
<script>
export default {
  emits: ['page'],
  props: {
    pages: Number,
    current: Number,
    start: {
      type: Number,
      default: 1,
    },
    end: {
      type: Number,
      default: 0,
    },
  },
  computed: {
    last: (c) => (c.end > 0 ? c.end : c.pages),
    parsed_pages: (c) => {
      const dt = [];
      for (let index = c.start; index <= c.last; index++) {
        dt.push({
          label: index,
          current: c.current === index,
        });
      }
      return dt;
    },
  },
};
</script>