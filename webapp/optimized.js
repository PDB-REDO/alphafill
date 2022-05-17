import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

import { Viewer } from "./molstar";

window.addEventListener('load', () => {
	const loadAlert = document.getElementById('load-alert');
	const errorAlert = document.getElementById('error-alert');
	const viewer = new Viewer(document.getElementById('app'));

	fetch(`v1/aff/${AF_ID}/optimized-with-stats/A,${ASYM_ID}`)
		.then(data => {
			loadAlert.style.display = 'none';
			return data.json();
		})
		.then(result => {
			viewer.loadStructureFromString(result.model)
				.then(() => {
					viewer.selectAsym(ASYM_ID);
				});

			const formatter = new Intl.NumberFormat('en-US', {
				minimumFractionDigits: 2,
				maximumFractionDigits: 2,
				});

			document.getElementById("combined-original").textContent = formatter.format(+result.lev.before.combined);
			document.getElementById("combined-optimized").textContent = formatter.format(+result.lev.after.combined);
			document.getElementById("protein-original").textContent = formatter.format(+result.lev.before.protein);
			document.getElementById("protein-optimized").textContent = formatter.format(+result.lev.after.protein);
			document.getElementById("ligand-original").textContent = formatter.format(+result.lev.before.ligand);
			document.getElementById("ligand-optimized").textContent = formatter.format(+result.lev.after.ligand);

			document.getElementById('model').classList.remove('invisible');
		})
		.catch(err => {
			loadAlert.style.display = 'none';
			errorAlert.style.display = '';
			alert(err);
		});
})
