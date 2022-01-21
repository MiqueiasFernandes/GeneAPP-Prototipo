<template>
  Single View genes carregados {{ genes.length }}
  <Button :disabled="!has_prev" @click="prev()">PREV</Button>
  <Button :disabled="!has_next" @click="next()">NEXT</Button>
  <ol v-for="file of files" :key="file.name">
    <li>{{ file.name }}</li>
  </ol>

  <div v-if="gene">
    gene ativo {{ gene.nome }} {{ gene.sequencia }} {{ gene.inicio }}
    {{ gene.fim }}
    <div id="canvas"></div>
  </div>
</template>

<script>
import Drawable from "../core/d3/Drawable";
import Bounds from "../core/d3/Bounds";
import DrawableGene from "../core/d3/DrawableGene";
import Histogram from "../core/d3/Histogram";

export default {
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
      genes: [],
      gene: null,
      drawable: null,
      files: [],
    };
  },
  mounted() {},
  methods: {
    set_data(genes, files) {
      console.log(genes);
      console.log(files);
      genes.forEach((g) => this.$data.genes.push(g));
      files.forEach((f) => this.$data.files.push(f));
      this.gene = this.genes[0];
      this.plotarGene();
    },
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
      const dg = new DrawableGene(
        this.drawable,
        this.gene,
        new Bounds(600, 300, 0, 60, { left: 82, right: 15 })
      );
      dg.draw();

      const dt = {}
this.files[0].data.map(x => ([parseInt(x[1]), parseInt(x[3])])).filter(a => a[0] > 0)
.forEach(x => dt[x[0]] =x[1] );
this.files[0].data.map(x => ([parseInt(x[2]), parseInt(x[3])])).filter(a => a[0] > 0)
.forEach(x => dt[x[0]] = (x[1] > 0 ? Math.log2(x[1]) : 0) );
console.log(dt)
      const hist = new Histogram(
        dg,
        dg.bounds
          .down(200)
          .withX(dg.bounds.x)
          .withWidth(dg.bounds.width)
          .withHeight(100),
          {
            samp1: { rep1:  dt},
            samp2: {rep1: {}}
          }
      );
      hist.draw();
    },
  },
};
</script>

<style>
</style>