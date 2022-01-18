<template>
  Single View genes carregados {{ genes.length }}
  <Button :disabled="!has_prev" @click="prev()">PREV</Button>
  <Button :disabled="!has_next" @click="next()">NEXT</Button>

  <div v-if="gene">
    gene ativo {{ gene.nome }}
    <div id="canvas"></div>
  </div>
</template>

<script>
import Drawable from "../core/d3/Drawable";
import Bounds from "../core/d3/Bounds";
import DrawableGene from "../core/d3/DrawableGene";

export default {
  props: {
    genes: {
      type: Array,
      default: () => [],
    },
  },
  computed: {
    has_prev() {
      return this.genes.indexOf(this.gene) > 0;
    },
    has_next() {
      return this.genes.length - this.genes.indexOf(this.gene) > 1;
    },
  },
  data() {
    return {
      gene: null,
      drawable: null,
    };
  },
  mounted() {
    this.gene = this.genes[0];
    this.plotarGene();
  },
  methods: {
    prev() {
      this.gene = this.genes[this.genes.indexOf(this.gene) - 1];
      this.plotarGene();
    },
    next() {
      this.gene = this.genes[this.genes.indexOf(this.gene) + 1];
      this.plotarGene();
    },

    clear() {
      if (this.drawable) this.drawable.clear("canvas");
    },

    plotarGene() {
      this.clear();
      this.drawable = new Drawable(
        null,
        "canvas",
        new Bounds(800, 600, 0, 0, {
          top: 20,
          bottom: 50,
          right: 200,
          left: 100,
        })
      );
      new DrawableGene(
        this.drawable,
        this.gene,
        new Bounds(600, 300, 0, 60, { left: 82, right: 15 })
      ).draw();
    },
  },
};
</script>

<style>
</style>