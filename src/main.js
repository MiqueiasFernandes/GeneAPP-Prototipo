import { createApp } from 'vue'
import bootstrap from './plugins/bootstrap'
import App from './App.vue'

const app = createApp(App).use(bootstrap)

app.config.globalProperties.$bootstrap_icons.then(() => {
    app.mount('#app')
    console.log('\n\
     ██████╗ ███████╗███╗   ██╗███████╗ █████╗ ██████╗ ██████╗ \n\
    ██╔════╝ ██╔════╝████╗  ██║██╔════╝██╔══██╗██╔══██╗██╔══██╗\n\
    ██║  ███╗█████╗  ██╔██╗ ██║█████╗  ███████║██████╔╝██████╔╝\n\
    ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝  ██╔══██║██╔═══╝ ██╔═══╝ \n\
    ╚██████╔╝███████╗██║ ╚████║███████╗██║  ██║██║     ██║    \n\
     ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝    \n\
    by bio@mikeias.net\n')
})
