import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

/* global hash */

var timer;

function update()
{
	fetch(`v1/aff/CS-${hash}/status`)
		.then(r => r.json())
		.then(r => {
			const statusContainer = document.getElementById('status-container');
			
			switch (r.status) {
				case "queued":
					statusContainer.classList.remove("unknown");
					statusContainer.classList.add("queued");
					break;
				case "running":
					statusContainer.classList.remove("unknown", "queued");
					statusContainer.classList.add("running");
					const progressBar = document.querySelector(".progress-bar");
					progressBar.style.width = (r.progress * 100) + "%";
					break;
				case "finished":
					statusContainer.classList.remove("unknown", "queued", "running");
					window.location = `model?id=CS-${hash}`;
					clearTimeout(timer);
				default:
					clearTimeout(timer);
			}
		});
}

window.addEventListener('load', () => {

	timer = setInterval(update, 1000);

});