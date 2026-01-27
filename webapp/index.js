import "core-js/stable";
import 'bootstrap';

window.addEventListener('load', () => {

	const uploadBtn = document.getElementById('upload-btn');
	if (uploadBtn) {
		uploadBtn.addEventListener('click', async (evt) => {
			evt.preventDefault();

			const customFile = document.getElementById('custom-file');

			if (customFile.files.length == 1) {

				const fd = new FormData();

				fd.append("structure", customFile.files[0]);

				const r = await fetch("v1/aff", {
					'Accept': 'application/json',
					'method': "POST",
					'body': fd
				});
				const data = await r.json();

				if (r.ok)
					window.location = `model?id=${data.id}`;
				else if (typeof (data.error) === "string")
					alert(data.error);
				else
					alert(`Failed to upload data: ${r.statusText}`);
			}
		});
	}

});