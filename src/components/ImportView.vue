<template>

  <div class="accordion accordion-flush" id="accordionFlushExample">
    <div class="accordion-item">
      <h2 class="accordion-header" id="flush-headingOne">
        <button
          class="accordion-button"
          type="button"
          data-bs-toggle="collapse"
          data-bs-target="#flush-collapseOne"
          aria-expanded="false"
          aria-controls="flush-collapseOne"
        >
         <Icon class="me-3 success" sm name='check2-circle' v-if="tabela_as.length > 0"></Icon> 
         <Badge class="me-3" v-if="tabela_as.length < 1">1</Badge>  
         <strong>Load Analise data</strong>
        </button>
      </h2>
      <div
        id="flush-collapseOne"
        class="accordion-collapse collapse show"
        aria-labelledby="flush-headingOne"
        data-bs-parent="#accordionFlushExample"
      >
        <div class="accordion-body">


 <div class="corpo">
<Display subtitle ico="hammer" class="secao">Experiment settings</Display>
<div class="row">
  <div class="col-md-6">
<div class="form-floating mb-3">
  <input  class="form-control" id="floatingInput" placeholder="Controle" v-model="lab_controle" :disabled="tabela_as.length > 0">
  <label for="floatingInput">Controle</label>
</div>
  </div>
  <div class="col-md-4">
    <label for="range1" class="form-label">Quantidade de replicas: {{replica_1}}</label>
<input type="range" class="form-range" min="1" max="9" id="range1" v-model="replica_1" :disabled="tabela_as.length > 0">
  </div>
  <div class="col-md-2">
    <label for="exampleColorInput" class="form-label">Color</label>
<input type="color" class="form-control form-control-color" id="exampleColorInput" v-model="cor_controle" title="Choose your color" :disabled="tabela_as.length > 0">
  </div>
</div>

<div class="row">
  <div class="col-md-6">
<div class="form-floating mb-3">
  <input  class="form-control" id="floatingInput" placeholder="Tratamento" v-model="lab_tratamento" :disabled="tabela_as.length > 0">
  <label for="floatingInput">Tratamento</label>
</div>
  </div>
  <div class="col-md-4">
    <label for="range2" class="form-label">Quantidade de replicas: {{replica_2}}</label>
<input type="range" class="form-range" min="1" max="9" id="range2" v-model="replica_2" :disabled="tabela_as.length > 0">
  </div>
  <div class="col-md-2">
    <label for="exampleColorInput" class="form-label">Color</label>
<input type="color" class="form-control form-control-color" id="exampleColorInput"  v-model="cor_tratamento" title="Choose your color" :disabled="tabela_as.length > 0">
  </div>
</div>
 </div>



 <div class="corpo">
<Display subtitle ico="table" class="secao">Open table with DAS genes list</Display>
<p>
  Abrir uma tabela csv (separado por comma) de genes com eventos AS significantes, contendo as seguintes colunas:
  <ol>
    <li>Nome do gene {{ genes_lab }}</li>
    <li>ΔPSI do evento</li>
    <li>Significancia estatistica (FDR)</li>
  </ol>
  <File
    @parsed="import_tabela_eventos($event)"
    txt
    :types="['csv']" :opts="false"
    :max="1" :disabled="tabela_as.length > 0"
  ></File>
Genes redundantes serao filtrados pelo maior abs(ΔPSI) significante.
</p>
 </div>



 <div class="corpo">
<Display subtitle class="secao" ico="body-text">Load sequence data</Display>
<ul class="list-group list-group-flush" v-if="tabela_as.length > 0">
  
  
  
  <li class="list-group-item d-flex justify-content-between align-items-start">
    <div class="ms-2 me-auto">
      <div class="fw-bold">Carregar o arquivo GFF3 dos genes</div>
      
      
     <div class="form-check form-switch">
  <input class="form-check-input" type="checkbox" role="switch" id="flexSwitchCheckChecked" v-model="usar_api_ncbi">
  <label class="form-check-label" for="flexSwitchCheckChecked">Importar pela API do NCBI</label>
</div>

<div v-if="usar_api_ncbi" class="row d-flex justify-content-between align-items-center">
<div class="col">
  <Button @click="carregar_gff_ncbi()" ico="cloud-arrow-down">Carregar</Button></div>
<div class="col">
  <ProgressBar
    ref="progress"
    label
    :min="0"
    :max="genes.length" style="width:20rem"
  ></ProgressBar></div>

</div>
<div v-if="!usar_api_ncbi">
  <File
    @parsed="import_tabela_gff($event)"
    txt
    :types="['gff', 'gff3']" :opts="false"
    :max="1" :size="1000000"
  ></File>
</div>


    </div>
    <span class="badge bg-success rounded-pill me-1" v-if="gff_genes > 0">{{gff_genes}} genes</span>
    <span class="badge bg-success rounded-pill me-1" v-if="gff_mrnas > 0">{{gff_mrnas}} mRNAS</span>
    <span class="badge bg-success rounded-pill me-1" v-if="gff_exons > 0">{{gff_exons}} exons</span>
    <span class="badge bg-danger rounded-pill me-1" v-for="o of gff_organisms" :key="o"><i>{{o}}</i></span>
  </li>


  <li class="list-group-item d-flex justify-content-between align-items-start">
    <div class="ms-2 me-auto">
      <div class="fw-bold">Fasta dos genes</div>
          
     <div class="form-check form-switch">
  <input class="form-check-input" type="checkbox" role="switch" id="flexSwitchCheckChecked2" v-model="usar_api_ncbi_fasta" :disabled="gff_exons < 1">
  <label class="form-check-label" for="flexSwitchCheckChecked2">Importar pela API do NCBI</label>
</div>
<div v-if="usar_api_ncbi_fasta" class="row d-flex justify-content-between align-items-center">
<div class="col">
  <Button @click="carregar_fasta_ncbi()" ico="cloud-arrow-down" :disabled="gff_exons < 1">Carregar</Button></div>
<div class="col">
  <ProgressBar
    ref="progress_fasta"
    label
    :min="0"
    :max="genes.length" style="width:20rem" :disabled="gff_exons < 1"
  ></ProgressBar></div>

</div>
<div v-if="!usar_api_ncbi_fasta">
  <File :disabled="gff_exons < 1"
    @parsed="import_fasta_genes($event)"
    txt
    :types="['fasta', 'fa', 'fna']" :opts="false"
    :max="genes.length" :size="1000000"
  ></File>
</div>
    </div>
    <span class="badge bg-success rounded-pill" v-if="genes_fasta.length > 0">{{genes_fasta.length}}</span>
  </li>
  
  
  
  <li class="list-group-item d-flex justify-content-between align-items-start">
    <div class="ms-2 me-auto">
      <div class="fw-bold">Fasta das protein</div>

    </div>
    <span class="badge bg-primary rounded-pill">14</span>
  </li>
</ul>




 </div>



 <div class="corpo">
<Display ico="lightning" subtitle class="secao">Generate missed data</Display>

<ul class="list-group list-group-flush" v-if="tabela_as.length > 0 && tabela_as.length === genes_fasta.length">


    <li class="list-group-item d-flex justify-content-between align-items-start">
    <div class="ms-2 me-auto">
      <div class="fw-bold">Sequencia das isoformas</div>
      <Button @click="extrair_iso()">obter</Button>
    </div>
    <span class="badge bg-primary rounded-pill">{{exs}}</span>
  </li>


    <li class="list-group-item d-flex justify-content-between align-items-start">
    <div class="ms-2 me-auto">
      <div class="fw-bold">Categoria dos eventos</div>
 <Button @click="classificar_evt()">classificar</Button>
    </div>
    <span class="badge bg-primary rounded-pill">{{as_evts.length}}</span>
  </li>


    <li class="list-group-item d-flex justify-content-between align-items-start">
    <div class="ms-2 me-auto">
      <div class="fw-bold">Sequencia das juncoes</div>
 <Button @click="sequencia_das_juncoes()">extrair</Button>
    </div>
    <span class="badge bg-primary rounded-pill">14</span>
  </li>



</ul>
 </div>








  <div v-if="carregado">
    <ol v-for="gene of genes" :key="gene.name">
      <li>
        <b>{{ gene.name }} {{ gene.status }}</b
        ><br />
        {{ gene.getParsed() }}
      </li>
    </ol>
  </div>


        </div>
      </div>
    </div>
    <div class="accordion-item">
      <h2 class="accordion-header" id="flush-headingTwo">
        <button
          class="accordion-button collapsed"
          type="button"
          data-bs-toggle="collapse"
          data-bs-target="#flush-collapseTwo"
          aria-expanded="false"
          aria-controls="flush-collapseTwo"
        >
          <Badge class="me-3">2</Badge>  <strong>Load transcriptomic data</strong>
        </button>
      </h2>
      <div
        id="flush-collapseTwo"
        class="accordion-collapse collapse"
        aria-labelledby="flush-headingTwo"
        data-bs-parent="#accordionFlushExample"
      >
        <div class="accordion-body">
          <ul>
            <li>bed</li>
            <li>bed juncao</li>
          </ul>
<p>
<span>Carregar bed da juncao</span>
  <File
    @parsed="import_bed_juncao($event)"
    txt
    :types="['bed']" :opts="false"
    :max="1" :size="1000000"
  ></File>   
</p>

  <File
    @files="ver($event)"
    @load="prg($event)"
    @parsed="read($event)"
    txt
    :types="['bed']"
    :max="10"
  ></File>

  <ProgressBar ref="prog" labelabs :min="0" :max="files.length"></ProgressBar>

        </div>
      </div>
    </div>
    <div class="accordion-item">
      <h2 class="accordion-header" id="flush-headingThree">
        <button
          class="accordion-button collapsed"
          type="button"
          data-bs-toggle="collapse"
          data-bs-target="#flush-collapseThree"
          aria-expanded="false"
          aria-controls="flush-collapseThree"
        >
         <Badge class="me-3">3</Badge> <strong>  Load extra data </strong>
        </button>
      </h2>
      <div
        id="flush-collapseThree"
        class="accordion-collapse collapse"
        aria-labelledby="flush-headingThree"
        data-bs-parent="#accordionFlushExample"
      >
        <div class="accordion-body">
          <ul>
            <li>tracks</li>
            <li>anotacao interpro (baixar | api)</li>
          </ul>
        </div>
      </div>
    </div>
  </div>

  <!-- <a href="coverage.bed">download BED de exemplo</a> -->




</template>

<script>
import NCBI from "../core/api/NCBI";
import AS from "../core/locus/AS";
// import Interpro from "../core/api/Interpro";
import Fasta from "../core/locus/Fasta";
// import Gene from '../core/locus/Gene';

export default {
  emits: ["pronto"],
  computed: {
    genes_lab(t) {
      return t.tabela_as.length > 0
        ? `(${[...new Set(t.tabela_as)].length} genes carregados)`
        : "";
    },
  },
  data() {
    return {
      carregado: false,

      replica_1: 3,
      replica_2: 3,
      lab_controle: null,
      lab_tratamento: null,
      cor_controle: "#B9D3EE",
      cor_tratamento: "#6941A4",
      usar_api_ncbi: false,
      gff_genes: 0,
      gff_mrnas: 0,
      gff_exons: 0,
      gff_organisms: [],
      usar_api_ncbi_fasta: false,
      genes_fasta: [],
      exs: 0,
      as_evts: [],
      bed_juncao: [],

      genes: [],
      files: [],
      tabela_as: [],
    };
  },
  mounted() {
    this.lab_controle = this.storage.get(`lab_controle`, this.lab_controle);
    this.lab_tratamento = this.storage.get(
      `lab_tratamento`,
      this.lab_tratamento
    );
    this.replica_1 = this.storage.get(`replica_1`, this.replica_1);
    this.replica_2 = this.storage.get(`replica_2`, this.replica_2);
    this.cor_controle = this.storage.get(`cor_controle`, this.cor_controle);
    this.cor_tratamento = this.storage.get(
      `cor_tratamento`,
      this.cor_tratamento
    );
    this.usar_api_ncbi = this.storage.get(`email`, "no").includes("@");
    this.usar_api_ncbi_fasta = this.storage.get(`email`, "no").includes("@");
  },
  methods: {
    import_tabela_eventos(tabela) {
      this.tabela_carregada = true;
      let data = tabela[0].txt
        .replaceAll('"', "")
        .replaceAll("'", "")
        .split("\n")
        .filter((l) => l.includes(","))
        .map((l) => l.split(","))
        .map((l) => [l[0], parseFloat(l[1]), parseFloat(l[2])])
        .slice(1);
      const nova = {};
      data.forEach((x) => {
        if (nova[x[0]]) {
          if (Math.abs(x[1]) > Math.abs(nova[x[0]][1])) {
            nova[x[0]] = x;
          }
        } else {
          nova[x[0]] = x;
        }
      });
      data = Object.values(nova);
      this.tabela_as = data;
      this.$toast(
        `${[...new Set(data)].length} genes carregados`,
        "OK!",
        "success",
        "check"
      );
      this.storage.set(`lab_controle`, this.lab_controle);
      this.storage.set(`lab_tratamento`, this.lab_tratamento);
      this.storage.set(`replica_1`, this.replica_1);
      this.storage.set(`replica_2`, this.replica_2);
      this.storage.set(`cor_controle`, this.cor_controle);
      this.storage.set(`cor_tratamento`, this.cor_tratamento);
    },
    carregar_gff_ncbi() {
      this.genes = NCBI.api().async_load(this.tabela_as.map((x) => x[0]));
      const ti = setInterval(() => {
        const qtd = this.genes.filter((g) => g.status === "loaded").length;
        this.$refs.progress.setValue(qtd);
        if (qtd === this.genes.length) {
          clearInterval(ti);
          this.gff_loaded();
        }
      }, 1000);
    },

    gff_loaded() {
      this.genes = NCBI.api().getGenes(this.genes);
      this.gff_genes = this.genes.length;
      this.gff_mrnas = this.genes
        .map((g) => g.getIsoformas().length)
        .reduce((a, b) => a + b);
      this.gff_exons = this.genes
        .map((g) => {
          if (g.getIsoformas().length < 1) console.log(g.nome);
          return g;
        })
        .map((g) =>
          g
            .getIsoformas()
            .map((i) => i.getExons().length)
            .reduce((a, b) => a + b)
        )
        .reduce((a, b) => a + b);
      this.gff_organisms = [...new Set(this.genes.map((g) => g.organism))];
    },
    import_tabela_gff(file) {
      // const gff = file[0].txt.split("\n").map(x => x.split("\t")).filter(l => !l[0].startsWith("#"));
      // const genes_interesse = this.tabela_as.map((x) => x[0]);
      // const genes_interesse_ncbi = genes_interesse.map((x) => 'gene-'+x);

      // const genes = gff
      // .filter(l => l[2] === 'gene')
      // .map(g => [g[8].split("ID=")[1].split(";")[0], g])
      // .map(g => [g[0], g[1], genes_interesse.includes(g[0]), genes_interesse_ncbi.includes(g[0])])
      // .filter(g => g[2] || g[3])

      // const mrnas = gff
      // .filter(l => l[2] === 'mRNA' && l[8].includes("Parent="))
      // .map(m => [m[8].split("Parent=")[1].split(";")[0], m[8].split("ID=")[1].split(";")[0], m])
      // .filter(m => genes_interesse.includes(m[0]) || genes_interesse_ncbi.includes(m[0]))

      // const mrnas_ids = mrnas.map(m => m[1])
      // const exons = gff
      // .filter(l => l[2] === 'exon' && l[8].includes("Parent="))
      // .map(e => [e[8].split("Parent=")[1].split(";")[0], e])
      // .filter(e => mrnas_ids.includes(e[0]))

      // const n_genes = Object.fromEntries(genes.map(g =>
      // [
      //   g[2] ? g[0] : g[0].remove('gene-'),
      //   new Gene(g[2] ? g[0] : g[0].remove('gene-'), parseInt(g[1][3]), parseInt(g[1][4])), g[1][4] === '+', g[1][0]
      // ]));
      console.log(`importar file ${file.name}`);
    },

    carregar_fasta_ncbi() {
      this.genes.forEach((g) =>
        NCBI.api()
          .get_locus(g.sequencia, g.inicio, g.fim, g.fita, g.nome)
          .then((f) => {
            if (g.size === f.size) {
              this.genes_fasta.push(f);
              this.$refs.progress_fasta.setValue(this.genes_fasta.length);
            } else {
              console.log(g);
              console.log(f);
            }
          })
      );
    },

    import_fasta_genes(files) {
      const targets = this.genes.map((g) => [
        `${g.sequencia}:${
          g.fita ? g.inicio + "-" + g.fim : "c" + g.fim + "-" + g.inicio
        }`,
        g,
      ]);
      files.forEach((file) => {
        Fasta.fromMultiFasta(file.txt)
          .map((x) => [
            targets.filter((t) => t[0] === x.id),
            this.genes.filter((g) => g.id === x.id),
            x,
          ])
          .filter((x) => x[0].length > 0 || x[1].length > 0)
          .forEach((x) => {
            if (x[0].length > 0) {
              if (x[0][0][1].size === x[2].size) {
                x[2].id = x[0][0][1].nome;
                this.genes_fasta.push((x[0][0][1] = x[2]));
              } else {
                console.log(x[0][0][1]);
                console.log(x[2]);
              }
            } else {
              if (x[1][0].size === x[2].size) {
                this.genes_fasta.push((x[1][0] = x[2]));
              } else {
                console.log(x[1][0]);
                console.log(x[2]);
              }
            }
          });
      });

      setTimeout(() => {
        if (this.genes_fasta.length !== this.genes) {
          const carr = this.genes_fasta.map((x) => x.id);
          const falt = this.genes
            .filter((g) => !carr.includes(g.nome))
            .map((g) => g.nome);
          this.$toast(
            `Faltam ${falt.length} genes agora ${
              falt.length < 10 ? falt.join(", ") : ""
            }`,
            "Fasta Carregado!",
            "success",
            "check"
          );
        }
      }, 3000);
    },

    extrair_iso() {
      this.genes.forEach((g) => {
        const fasta = this.genes_fasta.filter((f) => f.id === g.nome)[0];
        g.getIsoformas().forEach((i) => {
          i.getExons().forEach((e) => e.extrair(g, fasta.seq));
          this.exs++;
        });
      });
    },

    classificar_evt() {
      this.as_evts = this.genes.map((g) => new AS(g));
    },

    sequencia_das_juncoes() {
      /// juncoes de se
      this.as_evts.forEach((e) =>
        e.getSE().forEach((se) => {
          const gene = e.gene;
          const seq = this.genes_fasta.filter((f) => f.id === gene.nome)[0].seq;
          const exon_a1 = se.exon_canonico;
          const exon_a2 = se.exon;
          const exon_b1 = se.exons[0];
          const exon_b2 = se.exons[1];
          const ini_a1 = seq.indexOf(exon_a1.sequence);
          const ini_a2 = seq.indexOf(exon_a2.sequence);
          const ini_b1 = seq.indexOf(exon_b1.sequence);
          const ini_b2 = seq.indexOf(exon_b2.sequence);
          if (ini_a1 < 0 || ini_a2 < 0 || ini_b1 < 0 || ini_b2 < 0)
            console.log(gene);

          const menor = Math.min(
            ...[exon_a1.size, exon_a2.size, exon_b1.size, exon_b2.size, 100]
          );
          if (menor < 60) {
            /// se for menor pegar a seq no gene antes de ini_a1
            console.log(menor);
            console.log(gene);
          }

          const seqa1 = exon_a1.sequence.substr(
            exon_a1.sequence.length - menor
          );
          const seqa2 = exon_a2.sequence.substr(0, menor);
          const seqb1 = exon_b1.sequence.substr(
            exon_b1.sequence.length - menor
          );
          const seqb2 = exon_b2.sequence.substr(0, menor);
          se.canonic_seq = seqa1 + seqa2;
          se.alt_seq = seqb1 + seqb2;
        })
      );

      this.$emit("pronto", {
        genes: this.genes,
        as: this.as_evts,
        tabela: this.tabela_as,
        bed_juncao: this.bed_juncao
      });
    },

    import_bed_juncao(file) {
      const data = file[0].txt.split('\n').map(x => x.split('\t')).filter(x => x.length > 2);
      this.bed_juncao = data;
    },

    getGeneTable() {
      this.genes = NCBI.api().async_load(
        ///"LOC101258920,LOC101245278,LOC104644837,57502"
        [
          ...new Set(
            "LOC101252258,LOC101254260,LOC101247345,LOC101267654,LOC101246025,LOC101254143,LOC101254156,LOC101252360,LOC101268806,\
LOC101258125,LOC101261457,LOC101253816,LOC101259292,LOC104645871,LOC101261639,LOC101259814,LOC101252029,LOC104644837,\
LOC101263397,LOC101248986,LOC101255908,LOC101247787,LOC101268021,LOC101262589,LOC101258920,LOC101249485,SlPer1,LOC101249078,\
LOC101257882,LOC101264171,LOC101252559,LOC101266909,LOC101268591,LOC101245961,LOC101253757,544073,LOC101248300,LOC101256738,\
LOC101257984,LOC101257542,LOC101257100,LOC101256051,LOC101263341,LOC543863,LOC101246155,LOC101248222,LOC101252721,\
LOC101253077,LOC101267684,LOC101264788,LOC101246893,CYP85A1,LOC101257064,LOC101268345,LOC101259103,LOC101264545,\
LOC101245278,LOC101256296,LOC101245201,LOC101253775,LOC101254150,LOC101257084,LOC101246361,LOC101245882,LOC101251722,\
LOC101263627,LOC101267163,LOC101267838,LOC101251219,LOC101257931,LOC101247498,LOC101256522,LOC101244484,LOC101259888,\
LOC101252600,LOC101263245,LOC101263952,LOC101263285,LOC101259772,LOC101268326,LOC101248314,LOC101258184,LOC101262041,\
LOC101264109,LOC101256831,LOC101254006,LOC101260719,LOC101253917,LOC101246357,LOC101268268,LOC101265213,LOC101253473,\
LOC101244006,LOC101255628,LOC101260211,LOC101266673,LOC101245489,LOC101247954,LOC101261960,LOC101266993".split(
              ","
            )
          ),
        ]
      );

      const ti = setInterval(() => {
        const qtd = this.genes.filter((g) => g.status === "loaded").length;
        console.log(qtd);
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
              // const job = Interpro.api().post(
              //   fasta.id,
              //   fasta.seq,
              //   "bio@mikeias.net"
              // );
              // job.cbk = (x) => console.log(x);
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
.secao {
  margin-top: -2.6rem;
  background: white;
  width: 50%;
  padding-left: 1.5rem;
  margin-left: 2rem;
}
.corpo {
  padding: 1rem;
  border: 1px solid #adadad;
  border-radius: 0.5rem;
  margin-top: 2.5rem;
}
</style>