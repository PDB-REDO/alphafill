import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

import { Viewer } from "./molstar";

window.addEventListener('load', () => {

	const viewer = new Viewer(document.getElementById('app'), {
		layoutIsExpanded: false,
		layoutShowControls: false,
		viewportShowExpand: false,
		collapseLeftPanel: true,
		pixelScale: 1,
		pickScale: 0.25,
		pickPadding: 1,
	});

	viewer.loadStructureFromUrl(`/v1/aff/${AF_ID}`)
		.then(() => {
			const rows = document.querySelectorAll("tr.transplanted-row");
			[...rows].forEach(row => {
				row.addEventListener('click', () => {
					const asymID = row.getAttribute('data-asym-id');
					viewer.selectAsym(asymID);
				});
			});
		});
})
