<template>
  <div class="container mt-5">
    <GeneView></GeneView>
  </div>
  <Dialog global />
  <Toast global pos="'bottom'" />
</template>

<script>
import GeneView from "./components/GeneView.vue";

export default {
  name: "App",
  components: {
    GeneView,
  },
};
</script>