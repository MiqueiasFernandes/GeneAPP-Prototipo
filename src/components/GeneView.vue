<template>
  <Display><Icon :name="'pen'"></Icon>GeneView</Display>

  <div class="d-flex justify-content-evenly mt-3">
    <Button ico="arrow-left" disabled sm>Anterior</Button>

    <Button
      ico="info-circle"
      fill
      :color="'info'"
      outline
      sm
      v-if="gene"
      v-popover:bottom="{ title: 'Gene1', html: true, content: buildInfo }"
      >Informações</Button
    >
    <Button ico="rulers" @click="hide()" secondary sm>Régua</Button>
    <Button ico="cloud-arrow-down" fill @click="download()" secondary sm
      >Download</Button
    >
    <Button ico="share" fill @click="compartilhar()" success sm
      >Compartilhar</Button
    >
    <Button
      ico="chat-text"
      fill
      @click="comment()"
      secondary
      sm
      :alerta="comentario"
    >
      Comentário
    </Button>
    <Button
      :ico="star"
      @click="star = star == 'star' ? 'star-fill' : 'star'"
      sm
      warning
      outline
      >Boa</Button
    >
    <Button ico="arrow-right" sm>Próximo</Button>
  </div>
  <div id="canvas"></div>
</template>


<script>
import Isoforma from "../core/locus/Isoforma";
import Exon from "../core/locus/Exon";
import Locus from "../core/locus/Locus";
import CDS from "../core/locus/CDS";
import Gene from "../core/locus/Gene";
import Dominio from "../core/locus/Dominio";
import Drawable from "../core/d3/Drawable";
import Bounds from "../core/d3/Bounds";
import DrawableGene from "../core/d3/DrawableGene";

export default {
  setup() {},

  data() {
    return {
      regua: null,
      regua_status: true,
      star: "star",
      gene: null,
      comentario: "",
    };
  },

  computed: {
    buildInfo() {
      return [
        `<b>Cromossomo:</b> Chr4 ${this.gene.gene.fita ? "" : "anti-"}senso`,
        `<b>Coordenadas:</b> ${this.gene.gene.inicio}:${this.gene.gene.fim}`,
        ...this.gene.gene
          .getIsoformas()
          .map((i) => `<b>${i.nome}:</b> ${i.inicio}:${i.fim}`),
      ].join("<br>");
    },
  },

  mounted() {
    const real = () => {
      const fita = true;
      const gene = new Gene("LOC101245278", 63799797, 63804583, fita);
      const iso1 = new Isoforma("XM_019215223.2", 63799797, 63804583, fita);
      iso1.addExon(new Exon(null, 63799797, 63800877, fita));
      iso1.addExon(new Exon(null, 63803859, 63804583, fita));
      iso1.addCDS(new CDS(null, 63799996, 63800877, fita));
      iso1.addCDS(new CDS(null, 63803859, 63804467, fita));
      const iso2 = new Isoforma("XM_004245598.4", 63799813, 63802345, fita);
      iso2.addExon(new Exon(null, 63799813, 63800877, fita));
      iso2.addExon(new Exon(null, 63801138, 63802345, fita));
      iso2.addCDS(new CDS(null, 63799996, 63800877, fita));
      iso2.addCDS(new CDS(null, 63801138, 63801746, fita));
      const iso3 = new Isoforma("XM_010327128.3", 63801748, 63804583, fita);
      iso3.addExon(new Exon(null, 63801748, 63802340, fita));
      iso3.addExon(new Exon(null, 63802839, 63803616, fita));
      iso3.addExon(new Exon(null, 63803859, 63804583, fita));
      iso3.addCDS(new CDS(null, 63802906, 63803616, fita));
      iso3.addCDS(new CDS(null, 63803859, 63804467, fita));
      const iso4 = new Isoforma("XM_004245599.4", 63802499, 63804583, fita);
      iso4.addExon(new Exon(null, 63802499, 63803616, fita));
      iso4.addExon(new Exon(null, 63803859, 63804583, fita));
      iso4.addCDS(new CDS(null, 63802780, 63803616, fita));
      iso4.addCDS(new CDS(null, 63803859, 63804467, fita));
      gene.addIsoforma(iso1);
      gene.addIsoforma(iso2);
      gene.addIsoforma(iso3);
      gene.addIsoforma(iso4);
      return gene;
    };

    const mockGene = (nome = "Gene1", fita = true) => {
      const isoformas = {
        Isoforma1: [
          [200, 250],
          [280, 300],
          [360, 390],
          [480, 500],
        ],
        Isoforma2: [
          [270, 310],
          [350, 400],
          [450, 500],
        ],
        Isoforma3: [
          [240, 250],
          [260, 370],
          [420, 430],
        ],
      };

      const cds = {
        Isoforma1: [
          [280, 300],
          [360, 390],
        ],
        Isoforma2: [
          [380, 400],
          [450, 480],
        ],
        Isoforma3: [[330, 360]],
      };

      const isos = Object.keys(isoformas).map((i) => {
        const exons = isoformas[i].map(
          (e, idx) => new Exon("e" + idx, e[0], e[1], fita)
        );

        const isoforma = new Isoforma(
          i,
          Math.min(...exons.map((e) => e.inicio)),
          Math.max(...exons.map((e) => e.fim)),
          fita
        );
        exons.forEach((e) => isoforma.addExon(e));
        cds[i].forEach((c) => isoforma.addCDS(new CDS(null, c[0], c[1], fita)));
        return isoforma;
      });

      const gene = new Gene(
        nome,
        Math.min(...isos.map((e) => e.inicio)),
        Math.max(...isos.map((e) => e.fim)),
        fita
      );

      isos.forEach((iso) => gene.addIsoforma(iso));
      return gene;
    };

    const gene1 = mockGene();
    real();
    const gene2 = mockGene("Gene2", false);

    const cd1 = "Controle";
    const cd2 = "Tratamento";

    gene1.setExpressao(5, cd1);
    gene1.setExpressao(20, cd2);

    gene1.getIsoformas()[0].setExpressao(8, cd1);
    gene1.getIsoformas()[0].setExpressao(8, cd2);
    gene1.getIsoformas()[1].setExpressao(10, cd1);
    gene1.getIsoformas()[1].setExpressao(3, cd2);
    gene1.getIsoformas()[2].setExpressao(5, cd1);
    gene1.getIsoformas()[2].setExpressao(10, cd2);

    gene1.getIsoformas()[0].addDominio(new Dominio("dom1", 3, 40));
    gene1.getIsoformas()[1].addDominio(new Dominio("dom1", 5, 12));
    gene1.getIsoformas()[2].addDominio(new Dominio("dom1", 5, 9));

    gene1.getIsoformas()[0].addDominio(new Dominio("dom1", 450, 480));
    gene1.getIsoformas()[1].addDominio(new Dominio("dom1", 160, 350));
    gene1.getIsoformas()[2].addDominio(new Dominio("dom1", 30, 80));

    gene2.getIsoformas()[0].addDominio(new Dominio("dom1", 1, 5));
    gene2.getIsoformas()[1].addDominio(new Dominio("dom1", 9, 16));
    gene2.getIsoformas()[2].addDominio(new Dominio("dom1", 5, 8));

    const drawable = new Drawable(
      null,
      "canvas",
      new Bounds(800, 500, 0, 0, { top: 20, bottom: 50, right: 200, left: 100 })
    );

    const transposon_draw_strategy = (l) =>
      l.wave(l.bounds.x, l.bounds.y, l.bounds.width, "blue");
    gene1.addLocus(
      new Locus("transp", 220, 280, gene1.fita, "TRANSPOSON").set_draw_strategy(
        transposon_draw_strategy
      )
    );

    const snp_draw_strategy = (l) =>
      l.rect(l.bounds.x, l.bounds.y, l.bounds.width, 15);
    gene1.addLocus(
      new Locus("snp", 200, 200, gene1.fita, "SNP").set_draw_strategy(
        snp_draw_strategy
      )
    );
    gene1.addLocus(
      new Locus("snp", 300, 300, gene1.fita, "SNP").set_draw_strategy(
        snp_draw_strategy
      )
    );
    gene1.addLocus(
      new Locus("snp", 405, 405, gene1.fita, "SNP").set_draw_strategy(
        snp_draw_strategy
      )
    );
    gene1.addLocus(
      new Locus("snp", 408, 408, gene1.fita, "SNP").set_draw_strategy(
        snp_draw_strategy
      )
    );
    gene1.addLocus(
      new Locus("snp", 410, 410, gene1.fita, "SNP").set_draw_strategy(
        snp_draw_strategy
      )
    );
    gene1.addLocus(
      new Locus("snp", 500, 500, gene1.fita, "SNP").set_draw_strategy(
        snp_draw_strategy
      )
    );

    gene1.addLocus(
      new Locus("lncRNA", 220, 350, gene1.fita, "LNCRNA").set_draw_strategy(
        transposon_draw_strategy
      )
    );
    gene1.addLocus(
      new Locus("miRNA", 450, 460, gene1.fita, "MIRNA").set_draw_strategy(
        transposon_draw_strategy
      )
    );

    const juncao_draw_strategy = (l) => {
      l.text(
        l.bounds.mh(),
        l.bounds.y,
        `${l.locus.getExpressao(cd1)} / ${l.locus.getExpressao(cd2)}`,
        { hc: true, serif: true, fs: ".8em" }
      );
      l.curva({ x: l.bounds.x, y: l.bounds.y, x2: l.bounds.r });
      l.triangulo(l.bounds.r, l.bounds.y, 8);
    };

    gene1.addLocus(
      new Locus("SPLICING", 370, 430, gene1.fita, "SPLICING")
        .set_draw_strategy(juncao_draw_strategy)
        .setExpressao("500", cd1)
        .setExpressao("458.6", cd2)
    );

    this.$data.gene = new DrawableGene(
      drawable,
      gene1,
      new Bounds(600, 300, 0, 60, { left: 82, right: 15 })
    );
    this.$data.gene.draw();

    this.$data.regua = drawable.g(null, "regua");
    this.$data.regua.regua(
      gene1.inicio,
      gene1.fim,
      0,
      10,
      this.$data.gene.bounds.width,
      drawable.bounds.margin.left,
      null,
      true,
      600
    );
    this.hide();

    // new DrawableGene(
    //   drawable,
    //   gene2,
    //   new Bounds(400, 300, 20, 200, { right: 15, left: 82 })
    // ).draw();
  },

  methods: {
    hide() {
      if (this.$data.regua_status) {
        this.$data.regua.hide();
      } else {
        this.$data.regua.show();
      }
      this.$data.regua_status = !this.$data.regua_status;
    },

    download() {
      //this.$toast('Baixando...', 'title', 'success', 1000)
      this.$data.gene.download();
    },

    comment() {
      const salvar = (dt) => {
        this.$data.comentario = dt.comentario;
        return true;
      };
      const apagar = () => {
        this.$data.comentario = "";
        return true;
      };
      this.$dialog({
        sm: true,
        title: "Comentário",
        ico: "chat-text-fill",
        form: [
          {
            id: "comentario",
            hide_label: true,
            value: this.$data.comentario,
            placeholder: "digite aqui...",
          },
        ],
        actions: [
          {
            label: "Apagar",
            ico: "trash",
            color: "danger",
            click: (dt) => apagar(dt),
          },
          { label: "Salvar", ico: "check2", click: (dt) => salvar(dt) , auto: true}, 
        ],
      });
    },

    compartilhar(text) {
      text = text || window.location.href;
      const res = (title, txt, c) =>
        this.$toast(txt, title, c, "clipboard-check");
      navigator.clipboard.writeText(text).then(
        () => res("Copiado", "com sucesso para a área de transferencia!"),
        () => res("ERRO", "ao copiar para a área de transferencia.", "danger")
      );
    },
  },
};
</script>

<style scoped>
#canvas {
  margin-top: 1em;
  border: 1px solid gray;
  border-radius: 1.3em;
  text-align: center;
  background: rgb(248, 249, 250);
}
</style>