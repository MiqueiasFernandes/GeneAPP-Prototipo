<template>
  importt view

  <!-- <a href="coverage.bed">download BED de exemplo</a> -->

  <Button @click="getGeneTable()">carregar</Button>
  <File
    @files="ver($event)"
    @load="prg($event)"
    @parsed="read($event)"
    txt
    :types="['bed']"
    :max="10"
  ></File>

  <ProgressBar ref="prog" labelabs :min="0" :max="files.length"></ProgressBar>

  <ProgressBar
    v-if="genes.length > 0 && !carregado"
    ref="progress"
    label
    :min="0"
    :max="genes.length"
  ></ProgressBar>

  <div v-if="carregado">
    <ol v-for="gene of genes" :key="gene.name">
      <li>
        <b>{{ gene.name }} {{ gene.status }}</b
        ><br />
        {{ gene.getParsed() }}
      </li>
    </ol>
  </div>
</template>

<script>
import NCBI from "../core/api/NCBI";
import Interpro from "../core/api/Interpro";
import Fasta from "../core/locus/Fasta";

export default {
  emits: ["pronto"],
  data() {
    return {
      carregado: false,
      genes: [],
      files: [],
    };
  },
  mounted() {},
  methods: {
    getGeneTable() {
      this.genes = NCBI.api().async_load(
        "LOC101258920,LOC101245278,LOC104644837,57502".split(",")
      );

      const ti = setInterval(() => {
        const qtd = this.genes.filter((g) => g.status === "loaded").length;
        this.$refs.progress.setValue(qtd);
        if (qtd === this.genes.length || !this.$refs.progress) {
          clearInterval(ti);
          this.loaded();
        }
      }, 1000);
    },

    loaded() {
      this.carregado = true;
      const dt = [
        NCBI.api().getGenes(this.genes),
        this.files.map((f) => ({
          name: f.file.name,
          data: f.txt.split("\n").map((l) => l.split("\t")),
        })),
      ];
      dt[0].forEach((g) => {
        g.getIsoformas().forEach((i) =>
          NCBI.api()
            .get_sequence(i.getProteinName())
            .then((seq) => {
              const fasta = new Fasta(seq.data);
              console.log(fasta);
              const job = Interpro.api().post(
                fasta.id,
                fasta.seq,
                "bio@mikeias.net"
              );
              job.cbk = (x) => console.log(x);
            })
        );
      });
      console.log(dt);
      this.$emit("pronto", dt);
    },

    ver(files) {
      this.files = files;
    },

    prg(dt) {
      this.$refs.prog.setValue(dt[0]);
    },

    read(files) {
      this.files = files;
    },
  },
};
</script>

<style>
</style>