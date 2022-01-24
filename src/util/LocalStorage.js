export default class LocalStorage {

    static __INSTANCE = null;
    storage = null;
    prefix = "geneapp";

    static instance() {
        if (!LocalStorage.__INSTANCE)
            LocalStorage.__INSTANCE = new LocalStorage()
        LocalStorage.__INSTANCE.storage = window['localStorage'];
        return LocalStorage.__INSTANCE;
    }

    set(k, v) {
        this.storage.setItem(`${this.prefix}.${k.toLowerCase()}`, v);
    }

    get(k, or) {
        return this.storage.getItem(`${this.prefix}.${k.toLowerCase()}`) || or;
    }

    clear(k) {
        this.storage.removeItem(`${this.prefix}.${k.toLowerCase()}`);
    }

    storageAvailable(type = 'localStorage') {
        var storage;
        try {
            storage = window[type];
            var x = '__storage_test__';
            storage.setItem(x, x);
            storage.removeItem(x);
            return true;
        }
        catch (e) {
            return e instanceof DOMException && (
                // everything except Firefox
                e.code === 22 ||
                // Firefox
                e.code === 1014 ||
                // test name field too, because code might not be present
                // everything except Firefox
                e.name === 'QuotaExceededError' ||
                // Firefox
                e.name === 'NS_ERROR_DOM_QUOTA_REACHED') &&
                // acknowledge QuotaExceededError only if there's something already stored
                (storage && storage.length !== 0);
        }
    }

}


