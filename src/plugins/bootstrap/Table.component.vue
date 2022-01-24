<template>
  <table
    :class="{
      table: true,
      'table-sm': sm,
      'table-dark': dark,
      'table-striped': striped,
      'table-bordered': bordered,
      'thead-dark': head_dark,
      'table-hover': hover,
    }"
  >
    <thead>
      <tr>
        <th scope="col" v-for="head of heads" :key="head">{{ head }}</th>
      </tr>
    </thead>
    <tbody>
      <tr v-for="row of rows" :key="row[0]">
        <th scope="row">{{ row[0] }}</th>
        <td v-for="col of row.slice(1)" :key="col">{{ col }}</td>
      </tr>
    </tbody>
  </table>
</template>

<script>
export default {
  props: {
    sm: Boolean,
    dark: Boolean,
    striped: Boolean,
    bordered: Boolean,
    head_dark: Boolean,
    hover: Boolean,
  },
  data() {
    return {
      heads: [],
      rows: [],
    };
  },
  methods: {
      set_data(h, r) {
          this.heads = h
          this.rows = r
      }
  }
};
</script>

<style>
</style>