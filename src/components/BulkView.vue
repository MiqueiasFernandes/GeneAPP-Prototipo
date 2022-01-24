<template>
  <Tabs>
    <template #tab4-title> Tabela </template>
    <template #tab4>
      <Table ref="tabela" sm hover striped></Table>
      <ul>
        <li>
          tabela dos eventos
          chr,gene,iso,grupo,tipo,evento,psi,fdr,posa,posb,iso2,tpmgene1,tpmgene2,tpmiso1,tpmiso2,tpmiso21,
          tpmiso22,lengene,
          lenexon,lencds,lenintro,qtdmrna,qtdex,qtdin,emcds,cdsaltsize
        </li>
        <li>filtrar colunas</li>
        <li>filtrar por carac</li>
      </ul>
    </template>

    <template #tab5-title> samples </template>
    <template #tab5>
      <Button @click="drawViolin()">ver</Button>
      <div id="violin"></div>

      <ul>
        <li>
          violin do contexto genomico das features total, por amostra, por tipo
          evento heatmap evt/tipo/am
        </li>
      </ul>
    </template>

    <template #tab1-title>Junction reads coverage</template>
    <template #tab1>
      <Button @click="drawJuncao()">ver</Button>
      <div id="juncao"></div>

      <ul>
        <li>
          grafido de barras por amostra (juntar replicas) do tamanho da juncao
          amostrado empilhado por tipo basico de evento
        </li>
        <li>
          grafico de linha do tp basico e evt com a profundidade de reads
          down/up stream do sitio alternativo por amostra
        </li>
      </ul>
    </template>
    <template #tab2-title> Events overview </template>
    <template #tab2>
      <ul>
        <li>barras por tipo de evento</li>
        <li>volcano de evento, signifig por ∂psi</li>
      </ul>
    </template>

    <template #tab3-title> IsoStar </template>
    <template #tab3>
      <ul>
        <li>grafico de genes por amostra por iso com mais tpm</li>
      </ul>
    </template>
  </Tabs>
</template>

<script>
import Bounds from "../core/d3/Bounds";
import Drawable from "../core/d3/Drawable";
import Violin from "../core/d3/Violin";
import Line from "../core/d3/Line";
import File from "../util/File";

export default {
  data() {
    return {
      genes: [],
      as: [],
      tabela_heads: ["gene", "psi", "fdr"],
      tabela_rows: [],
      bed_juncao: [],
    };
  },
  mounted() {},
  methods: {
    set_data(dx) {
      this.genes = dx.genes;
      this.tabela = dx.tabela;
      this.as = dx.as;
      this.bed_juncao = dx.bed_juncao;

      this.$refs.tabela.set_data(["gene", "psi", "fdr"], this.tabela);
    },

    drawViolin() {
      const drawable = new Drawable(
        null,
        "violin",
        new Bounds(400, 300, 0, 0, {
          top: 20,
          bottom: 50,
          right: 100,
          left: 100,
        })
      );

      const lines = [];

      this.as.map((e) =>
        e.getSE().forEach((se) => {
          lines.push(
            `>S${lines.length} ${se.canonic_seq.length}\n${se.canonic_seq}`
          );
          lines.push(`>S${lines.length} ${se.alt_seq.length}\n${se.alt_seq}`);
        })
      );

      File.download("seqs.fa", lines.join("\n"));

      let dt = this.as
        .map((e) =>
          e
            .getSE()
            .map((x) => x.canonic_seq.length + "," + x.alt_seq.length)
            .join(",")
        )
        .join(",")
        .split(",");

      new Violin(drawable, drawable.bounds, {
        tamanho_exons_se: dt.map((x) => parseInt(x)),
      }).plot();
    },

    drawJuncao() {
      const drawable = new Drawable(
        null,
        "juncao",
        new Bounds(400, 300, 0, 0, {
          top: 20,
          bottom: 50,
          right: 100,
          left: 100,
        })
      );
      var data = {};
      this.bed_juncao.forEach((x) =>
        !data[x[0]]
          ? (data[x[0]] = [
              [(parseInt(x[1]) + parseInt(x[2])) / 2, parseInt(x[3])],
            ])
          : data[x[0]].push([
              (parseInt(x[1]) + parseInt(x[2])) / 2,
              parseInt(x[3]),
            ])
      );
      new Line(drawable, drawable.bounds, data).plot();
    },
  },
};
</script>

<style>
</style>