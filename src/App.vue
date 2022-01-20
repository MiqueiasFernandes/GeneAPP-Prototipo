<template>
  <div class="container">
    <header
      class="d-flex flex-wrap justify-content-center py-3 mb-4 border-bottom"
    >
      <a
        href="/"
        class="
          d-flex
          align-items-center
          mb-3 mb-md-0
          me-md-auto
          text-dark text-decoration-none
        "
      >
        <img src="./assets/img/logo.png" width="40" />
        <span class="fs-4">GeneAPP</span>
      </a>

      <ul class="nav nav-pills">
        <li class="nav-item">
          <a
            href="#"
            :class="{
              'nav-link': true,
              active: state < 2,
              disabled: state > 1,
            }"
            >Nova analise</a
          >
        </li>
        <li class="nav-item">
          <a
            href="#"
            @click="ir_pagina()"
            :class="{
              'nav-link': true,
              active: state == 2,
              disabled: !pronto || state > 2
            }"
            >Single View</a
          >
        </li>
      </ul>
    </header>
  </div>
  <div class="container mt-5">
    <ImportView @pronto="step1ok($event)" v-if="state == 1" ></ImportView>
    <!-- bulk view -->
    <SingleView v-if="state == 2" ref="single"></SingleView>
    <!-- download data -->
    <!-- sobre -->
  </div>
  <Dialog global />
  <Toast global pos="'bottom'" />
</template>

<script>
// import GeneView from "./components/GeneView.vue";
import ImportView from "./components/ImportView.vue";
import SingleView from './components/SingleView.vue';
require("./assets/img/logo.png");

export default {
  name: "App",
  components: {
    // GeneView,
    ImportView,
    SingleView
  },

  data() {
    return {
      genes: null,
      carregado: false,
      state: 1,
      pronto: false
    };
  },

  methods: {
    step1ok(data) {
      this.pronto = true;
      this.genes = data[0];
      this.files = data[1];
    },
    ir_pagina() {
      this.state++;
      setTimeout(() => {
        this.$refs.single.set_data(this.genes, this.files);
      }, 3000);
      
    },
  },
};
</script>