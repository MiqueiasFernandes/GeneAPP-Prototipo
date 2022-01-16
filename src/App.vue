<template>
 <div class="container">
    <header class="d-flex flex-wrap justify-content-center py-3 mb-4 border-bottom">
      <a href="/" class="d-flex align-items-center mb-3 mb-md-0 me-md-auto text-dark text-decoration-none">
        <img src="./assets/img/logo.png" width="40" >
        <span class="fs-4">GeneAPP</span>
      </a>

      <ul class="nav nav-pills">
        <li class="nav-item"><a href="#" class="nav-link active" aria-current="page">Novo</a></li>
        <li class="nav-item"><a href="#" class="nav-link">Features</a></li>
        <li class="nav-item"><a href="#" class="nav-link">Pricing</a></li>
        <li class="nav-item"><a href="#" class="nav-link">FAQs</a></li>
        <li class="nav-item"><a href="#" class="nav-link">About</a></li>
      </ul>
    </header>
  </div>
  <div class="container mt-5">
    <GeneView></GeneView>
  </div>
  <Dialog global />
  <Toast global pos="'bottom'" />
</template>

<script>
import GeneView from "./components/GeneView.vue";
require('./assets/img/logo.png');

export default {
  name: "App",
  components: {
    GeneView,
  },
};
</script>