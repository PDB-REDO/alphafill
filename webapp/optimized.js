import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

import { Viewer } from "./molstar";

window.addEventListener('load', () => {
	const loadAlert = document.getElementById('load-alert');
	const errorAlert = document.getElementById('error-alert');
	const viewer = new Viewer(document.getElementById('app'));

	fetch(`v1/aff/${AF_ID}/optimized/A,${ASYM_ID}/25`)
		.then(data => {
			loadAlert.style.display = 'none';
			return data.text();
		})
		.then(model => {
			viewer.loadStructureFromString(model);
			document.getElementById('model').classList.remove('invisible');
		})
		.catch(err => {
			loadAlert.style.display = '';
			alert(err);
		});
})
