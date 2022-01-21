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
        <li class="nav-item mx-2">
          <a
            href="#"
            :class="{
              'menu': true,
              'nav-link': true,
              active: state < 2,
              disabled: state > 1,
            }"
            ><Icon sm class="me-2" name='magic'></Icon>Import Data</a
          >
        </li>
        <li class="nav-item mx-2"><a href="#" class="nav-link"><Icon sm class="me-2" name='minecart-loaded'></Icon>
        Bulk View</a></li>
        
        <li class="nav-item mx-2">
          <a
            href="#"
            @click="ir_pagina()"
            :class="{
              'nav-link': true,
              active: state == 2,
              disabled: !pronto || state > 2
            }"
            ><Icon sm class="me-2" name='funnel'></Icon>Single View</a
          >
        </li>
        <li class="nav-item mx-2"><a href="#" class="nav-link"><Icon sm class="me-2" name='cloud-arrow-down'></Icon>Export Results</a></li>
        <li class="nav-item mx-2"><a href="#" class="nav-link"><Icon sm class="me-2" name='life-preserver'></Icon>About</a></li>
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
require("./assets/img/logo2.png");

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
<style scoped>
.menu {
  align-items: center;
    display: flex
}
</style>