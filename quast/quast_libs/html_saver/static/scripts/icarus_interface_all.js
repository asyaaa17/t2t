function setupInterface() {
    console.log("[setupInterface] Вызвана");

    if (typeof items === 'undefined') {
        return;
    }

    var checkboxes = document.getElementsByName('misassemblies_select');
    for (var i = 0; i < checkboxes.length; i++) {
        checkboxes[i].addEventListener('change', function () {
            showMisassemblies();
        });
    }

    if (checkboxes.length > 1) {
        checkboxes[checkboxes.length - 1].checked = false;
        showMisassemblies();
    }
}

window.setupInterface = setupInterface;