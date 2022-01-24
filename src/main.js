import { createApp } from 'vue'
import bootstrap from './plugins/bootstrap'
import App from './App.vue'


// Import the functions you need from the SDKs you need
import { initializeApp } from "firebase/app";
import { getAnalytics } from "firebase/analytics";
import LocalStorage from './util/LocalStorage';
// TODO: Add SDKs for Firebase products that you want to use
// https://firebase.google.com/docs/web/setup#available-libraries

// Your web app's Firebase configuration
// For Firebase JS SDK v7.20.0 and later, measurementId is optional
const firebaseConfig = {
  apiKey: "AIzaSyAaiFjUn293ACvOcpl0x3JCG2K_wk3-60I",
  authDomain: "geneapp-a6ab0.firebaseapp.com",
  projectId: "geneapp-a6ab0",
  storageBucket: "geneapp-a6ab0.appspot.com",
  messagingSenderId: "392143267162",
  appId: "1:392143267162:web:c65b3c7649e21531b8bf63",
  measurementId: "G-NDEHLDH3HS"
};

// // Initialize Firebase
const firebase_app = initializeApp(firebaseConfig);
// const analytics = 
getAnalytics(firebase_app);



const app = createApp(App).use(bootstrap);
app.config.globalProperties.storage = LocalStorage.instance();

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
