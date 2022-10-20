import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

/* global hash */

var timer;

function update() {
	var resultOK = true;
	fetch(`v1/aff/CS-${hash}/status`)
		.then(r => {
			resultOK = r.ok;
			return r.json();
		})
		.then(r => {
			const statusContainer = document.getElementById('status-container');

			statusContainer.classList.remove("unknown", "queued", "running", "error");

			if (resultOK == false) {
				statusContainer.classList.add("error");
				document.getElementById("errmsg").textContent = typeof (r.error) === "string" ? r.error : "Unknown error";
			}
			else {
				switch (r.status) {
					case "queued":
						statusContainer.classList.add("queued");
						break;
					case "running":
						statusContainer.classList.add("running");
						const progressBar = document.querySelector(".progress-bar");
						progressBar.style.width = (r.progress * 100) + "%";
						break;
					case "finished":
						window.location = `model?id=CS-${hash}`;
						clearTimeout(timer);
						break;
					case "unknown":
						statusContainer.classList.add("unknown");
						clearTimeout(timer);
						break;
					default:
						statusContainer.classList.add("error");
						document.getElementById("errmsg").textContent = typeof (r.message) === "string" ? r.message : "Unknown error";
						clearTimeout(timer);
				}

			}
		});
}

window.addEventListener('load', () => {

	timer = setInterval(update, 1000);

});