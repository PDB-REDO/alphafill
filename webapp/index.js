import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

// function readMyFile(file) {
//     return new Promise((resolve) => {
//         const fr = new FileReader();
//         fr.onload = e => resolve(e.target.result);
//         fr.readAsText(file);
//     })
// }

window.addEventListener('load', () => {

	const uploadBtn = document.getElementById('upload-btn');
	if (uploadBtn) {
		uploadBtn.addEventListener('click', () => {
			
			const customFile = document.getElementById('custom-file');
			const files = customFile.files;

			if (files.length == 1) {

				const fd = new FormData();

				fd.append("structure", files[0]);

				fetch("v1/aff", {
					'Accept': 'application/json',
					'method': "POST",
					'body': fd
				}).then(r => {
					if (r.ok)
						return r.json()
					throw "Failed to upload file";
				}).then(r => {
					window.location = `model?id=${r.id}`;
				}).catch(e => {
					console.log(e);
				});
			}
		 });
	}

});