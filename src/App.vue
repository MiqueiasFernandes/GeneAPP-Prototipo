<template>
  <div class="container">
    <header
      class="
        d-flex
        nav
        mainnav
        flex-wrap
        justify-content-center
        py-3
        mb-4
        border-bottom
        fixed-top
        shadow
        navbar-light
        px-5
      "
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
        <span class="fs-4 me-1">GeneAPP</span>
        <Badge color="info" round sm>1</Badge>
      </a>

      <ul class="nav nav-pills">
        <li class="nav-item mx-2">
          <a
            href="#"
            :class="{
              menu: true,
              'nav-link': true,
              active: state < 2,
              disabled: state > 1,
            }"
            ><Icon sm class="me-2" name="magic"></Icon>Import Data</a
          >
        </li>
        <li class="nav-item mx-2">
          <a
            href="#"
            :class="{
              'nav-link': true,
              active: state == 2,
              disabled: !pronto,
            }"
            @click="bulk()"
            ><Icon sm class="me-2" name="minecart-loaded"></Icon> Bulk View</a
          >
        </li>

        <li class="nav-item mx-2">
          <a
            href="#"
            @click="ir_pagina()"
            :class="{
              'nav-link': true,
              active: state == 3,
              disabled: !pronto,
            }"
            ><Icon sm class="me-2" name="funnel"></Icon>Single View</a
          >
        </li>
        <li class="nav-item mx-2">
          <a
            href="#"
            :class="{
              'nav-link': true,
              active: state == 4,
              disabled: !pronto,
            }"
            ><Icon sm class="me-2" name="cloud-arrow-down"></Icon>Export
            Results</a
          >
        </li>
        <!-- <li class="nav-item"><a href="#" class="nav-link"><Icon sm class="me-2" name='life-preserver'></Icon>About</a></li> -->

        <li class="nav-item mx-2 dropdown">
          <a
            class="nav-link dropdown-toggle"
            href="#"
            role="button"
            data-bs-toggle="dropdown"
            aria-expanded="false"
          >
            <Icon sm class="me-2" name="person"></Icon> More
          </a>
          <ul
            class="dropdown-menu"
            aria-labelledby="navbarDarkDropdownMenuLink"
          >
            <li>
              <a class="dropdown-item" href="#"
                ><Icon sm class="me-2" name="people"></Icon>About</a
              >
            </li>
            <li>
              <a class="dropdown-item" href="#"
                ><Icon sm class="me-2" name="life-preserver"></Icon>Help</a
              >
            </li>
            <li>
              <a
                class="dropdown-item"
                href="#"
                @click="(settings = !settings) && (outro = true)"
                ><Icon sm class="me-2" name="gear"></Icon>Settings</a
              >
            </li>
          </ul>
        </li>
      </ul>
    </header>
  </div>
  <div class="container content">
    <div class="card shadow-sm">
      <div class="card-body">
        <ImportView
          @pronto="step1ok($event)"
          v-if="state == 1 && !outro"
        ></ImportView>
        <BulkView ref="bulk" v-if="state == 2 && !outro"></BulkView>
        <SingleView v-if="state == 3 && !outro" ref="single"></SingleView>
        <!-- download data -->
        <!-- sobre -->

        <SettingsView v-if="settings"></SettingsView>
      </div>
    </div>
  </div>
  <Dialog global />
  <Toast global pos="'bottom'" />
</template>

<script>
// import GeneView from "./components/GeneView.vue";
import ImportView from "./components/ImportView.vue";
import BulkView from "./components/BulkView.vue";
import SingleView from "./components/SingleView.vue";
import SettingsView from "./components/SettingsView.vue";
require("./assets/img/logo2.png");

export default {
  name: "App",
  components: {
    // GeneView,
    ImportView,
    BulkView,
    SingleView,
    SettingsView,
  },

  data() {
    return {
      genes: null,
      carregado: false,
      state: 1,
      pronto: false,
      outro: false,
      settings: false,

      dx: null
    };
  },

  mounted() {
    const has_storage = this.storage.storageAvailable();

    if (!has_storage || this.storage.get("cookies", "no") !== "ok") {
      this.$toast(
        "This site use cookies",
        "Alert!",
        "secondary",
        "cone-striped",
        999999,
        () => has_storage && this.storage.set("cookies", "ok")
      );
    }

    if (!has_storage) {
      this.$toast(
        "Seu browser é incompativel com armazenamento local de dados.",
        "ERRO",
        "danger",
        "exclamation-triangle",
        30
      );
    } else {
      if (!this.storage.get("email")) {
        this.$toast(
          "Configure your user mail to use external APIs.",
          "Welcome!",
          "secondary",
          "flag-fill",
          30
        );
      }
    }
  },

  methods: {
    step1ok(data) {
      this.pronto = true;
      // this.genes = data[0];
      // this.files = data[1];
      this.dx = data;
    },
    bulk() {
      this.state = 2;
      setTimeout(() => {
        this.$refs.bulk.set_data(this.dx);
      }, 300);
    },
    ir_pagina() {
      this.state++;
      setTimeout(() => {
        this.$refs.single.set_data(this.dx);
      }, 3000);
    },
  },
};
</script>
<style scoped>
.mainnav {
  background: rgb(238 245 245 / 98%);
}
.menu {
  align-items: center;
  display: flex;
}
.content {
  margin-top: 6rem;
}
</style>