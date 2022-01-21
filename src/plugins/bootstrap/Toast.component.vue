<template>
  <div :class="['toasts', pos === 'top' ? 'onTop' : 'onBottom']">
    <div
      v-for="toast in toasts"
      :key="toast.id"
      :ref="toast.id"
      :id="toast.id"
      class="toast config d-flex align-items-center border-0"
      :class="[
        toast.color === 'light' ? '' : 'text-white',
        `bg-${toast.color}`,
      ]"
      role="alert"
      aria-live="assertive"
      aria-atomic="true"
    >
      <Icon v-if="toast.ico" :name="toast.ico" class="ms-2"></Icon>
      <div class="toast-body content">
        <strong class="mr-2 me-2" v-if="toast.title">{{ toast.title }}</strong>
        <span>{{ toast.text }}</span>
      </div>
      <button
        type="button"
        :class="toast.color === 'light' ? '' : 'btn-close-white'"
        class="btn-close ml-auto mr-2 fix-multiline"
        data-dismiss="toast"
        aria-label="Close"
        @click="toast.instance.hide()"
      ></button>
    </div>
  </div>
</template>
<script>
export default {
  props: {
    global: Boolean,
    pos: { type: String, default: "top" },
  },
  data: () => ({ toasts: [], index: 1 }),
  mounted() {
    if (this.global) {
      this.$toast_handler(this);
    }
  },
  updated() {
    this.toasts.forEach((toast) => {
      const el = this.$refs[toast.id];
      if (!toast.instance && el) {
        toast.instance = new this.$bootstrap.Toast(this.$refs[toast.id], {
          autohide: toast.delay > 0,
          delay: toast.delay * 1000,
        });
        toast.instance.show();
        el.addEventListener("shown.bs.toast", () => el.classList.add("shadow"));
        el.addEventListener("hidden.bs.toast", (el) => {
          const index = this.toasts.findIndex((t) => t.id === el.target.id);
          const toast_ = this.toasts[index];
          toast_.instance.hide();
          this.toasts.splice(index, 1);
        });
      }
    });
  },
  methods: {
    notify(text, title, color = "secondary", ico, delay = 5) {
      this.toasts.push({
        text,
        title,
        color,
        delay,
        ico,
        id: `toast-${this.index++}`,
      });
    },
  },
};
</script>
<style scoped>
.toasts {
  z-index: 100;
  position: fixed;
  right: 8rem;
}
.shadow {
  box-shadow: rgba(0, 0, 0, 0.25) 0px 54px 55px,
    rgba(0, 0, 0, 0.12) 0px -12px 30px, rgba(0, 0, 0, 0.12) 0px 4px 6px,
    rgba(0, 0, 0, 0.17) 0px 12px 13px, rgba(0, 0, 0, 0.09) 0px -3px 5px !important;
}
.onTop {
  top: 6rem;
}
.onBottom {
  bottom: 2rem;
}
.config {
  min-width: 500px;
  justify-content: space-between;
  padding: .5rem 1rem;
}

.content {
  width: 90%;
  display: flex;
  align-items: center;
}

.fix-multiline {
  padding: 0.4rem;
}
</style>