import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

window.addEventListener('load', () => {

	const uploadBtn = document.getElementById('upload-btn');
	if (uploadBtn) {
		uploadBtn.addEventListener('click', (evt) => {
			evt.preventDefault();

			const customFile = document.getElementById('custom-file');

			if (customFile.files.length == 1) {

				const fd = new FormData();

				fd.append("structure", customFile.files[0]);

				var resultOK = false;

				fetch("v1/aff", {
					'Accept': 'application/json',
					'method': "POST",
					'body': fd
				}).then(r => {
					resultOK = r.ok;
					return r.json()
				}).then(r => {
					if (resultOK) {
						window.location = `model?id=${r.id}`;
					}
					else if (typeof(r.error) === "string") {
						alert(r.error);
					}
					else {
						throw "Failed to upload file";
					}
				}).catch(e => {
					alert(e);
				});
			}
		});
	}

});